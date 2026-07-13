// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// -----------------------------------------------------------------------
// Coordinate descent update for one (i,j) element of D
// Mirrors CoordinateDescentUpdate() in quic.R
// -----------------------------------------------------------------------
static void coord_update(int i, int j, int J,
                         mat& D, const mat& Psi, const mat& Rho,
                         const mat& W, mat& U,
                         double& normD, double& diffD, const mat& S)
{
  double a = W(i,j) * W(i,j);
  if (i != j) a += W(i,i) * W(j,j);
  double ainv = 1.0 / a;

  // b = S[i,j] - W[i,j] + sum(W[i,] * U[,j])
  double b = S(i,j) - W(i,j) + dot(W.row(i), U.col(j).t());

  double c  = Psi(i,j) + D(i,j);
  double l  = Rho(i,j) * ainv;
  double f  = b * ainv;

  normD -= std::abs(D(i,j));

  double mu;
  if (c > f) {
    mu = -f - l;
    if (c + mu < 0.0) { mu = -c; D(i,j) = -Psi(i,j); }
    else               { D(i,j) += mu; }
  } else {
    mu = -f + l;
    if (mu + c > 0.0) { mu = -c; D(i,j) = -Psi(i,j); }
    else               { D(i,j) += mu; }
  }

  diffD += std::abs(mu);
  normD += std::abs(D(i,j));
  D(j,i) = D(i,j);

  if (mu != 0.0) {
    U.row(i) += mu * W.row(j);
    if (i != j) U.row(j) += mu * W.row(i);
  }
}

// -----------------------------------------------------------------------
// Diagonal Newton initialisation step
// Mirrors DiagNewton() in quic.R
// -----------------------------------------------------------------------
static double diag_newton(int J, const mat& S, const mat& Rho,
                           const mat& Psi, const mat& W, mat& D)
{
  // off-diagonal
  for (int i = 1; i < J; ++i) {
    for (int j = 0; j < i; ++j) {
      double a    = W(i,i) * W(j,j);
      double inva = 1.0 / a;
      double b    = S(i,j);
      double l    = Rho(i,j) * inva;
      double f    = b * inva;
      double mu;
      if (f < 0.0) {
        mu = -f - l;
        if (mu < 0.0) { mu = 0.0; D(i,j) = -Psi(i,j); }
        else           { D(i,j) += mu; }
      } else {
        mu = -f + l;
        if (mu > 0.0) { mu = 0.0; D(i,j) = -Psi(i,j); }
        else           { D(i,j) += mu; }
      }
    }
  }

  double logdet = 0.0, l1normX = 0.0, trSX = 0.0;

  for (int i = 0; i < J; ++i) {
    logdet  += std::log(Psi(i,i));
    l1normX += std::abs(Psi(i,i)) * Rho(i,i);
    trSX    += Psi(i,i) * S(i,i);

    double a    = W(i,i) * W(i,i);
    double ainv = 1.0 / a;
    double b    = S(i,i) - W(i,i);
    double l    = Rho(i,i) * ainv;
    double c    = Psi(i,i);
    double f    = b * ainv;
    double mu;
    if (c > f) {
      mu = -f - l;
      if (c + mu < 0.0) { D(i,i) = -Psi(i,i); continue; }
    } else {
      mu = -f + l;
      if (c + mu > 0.0) { D(i,i) = -Psi(i,i); continue; }
    }
    D(i,i) += mu;
  }

  D = D + D.t();
  D.diag() /= 2.0;

  return -logdet + trSX + l1normX;
}

// -----------------------------------------------------------------------
// Compute l1normX and trSX from Psi, S, Rho (lower triangle loop)
// -----------------------------------------------------------------------
static void compute_norms(int J, const mat& Psi, const mat& S, const mat& Rho,
                           double& l1normX, double& trSX)
{
  l1normX = 0.0; trSX = 0.0;
  for (int j = 0; j < J; ++j) {
    l1normX += Rho(j,j) * std::abs(Psi(j,j));
    trSX    += Psi(j,j) * S(j,j);
    for (int l = 0; l < j; ++l) {
      l1normX += 2.0 * Rho(j,l) * std::abs(Psi(j,l));
      trSX    += 2.0 * Psi(j,l) * S(j,l);
    }
  }
}

// -----------------------------------------------------------------------
// Main QUIC function — mirrors my_quic() in quic.R
// Returns a J x J x 2 array: [,,0] = Psi, [,,1] = W
// -----------------------------------------------------------------------
// [[Rcpp::export]]
arma::cube my_quic_cpp(int J, arma::mat S, arma::mat Rho,
                       double eps, int max_iter_quic)
{
  const double EPS_INNER  = 1e-8;
  const double cdSweepTol = 0.05;
  const int    max_lineiter = 20;
  const double sigma = 0.001;

  mat Psi = eye<mat>(J, J);
  mat W   = eye<mat>(J, J);
  mat D   = zeros<mat>(J, J);
  mat U   = zeros<mat>(J, J);

  // active set stored as pairs
  std::vector<int> as_row, as_col;
  as_row.reserve(J*(J+1)/2);
  as_col.reserve(J*(J+1)/2);

  double l1normX, trSX;
  compute_norms(J, Psi, S, Rho, l1normX, trSX);

  double logdetX = 0.0;
  double fX = 1e15, fX1 = 1e15, fXprev = 1e15;

  bool is_identity = true; // Psi starts as identity

  for (int newton = 0; newton < max_iter_quic; ++newton) {

    double normD = 0.0, diffD = 0.0, subgrad = 1e15;

    if (newton == 0 && is_identity) {
      D.zeros();
      fX = diag_newton(J, S, Rho, Psi, W, D);
    } else {
      int numActive = 0;
      U.zeros(); D.zeros(); subgrad = 0.0;
      as_row.clear(); as_col.clear();

      for (int j = 0; j < J; ++j) {
        for (int jj = 0; jj <= j; ++jj) {
          double g = S(j,jj) - W(j,jj);
          if (Psi(j,jj) != 0.0 || std::abs(g) > Rho(j,jj)) {
            as_row.push_back(j);
            as_col.push_back(jj);
            ++numActive;
            if      (Psi(j,jj) > 0.0) g += Rho(j,jj);
            else if (Psi(j,jj) < 0.0) g -= Rho(j,jj);
            else                       g  = std::abs(g) - Rho(j,jj);
            subgrad += std::abs(g);
          }
        }
      }

      int maxSweeps = (1 + newton) / 3;
      if (maxSweeps < 1) maxSweeps = 1;

      for (int sweep = 0; sweep < maxSweeps; ++sweep) {
        diffD = 0.0;
        // random shuffle of active set
        for (int i = numActive - 1; i > 0; --i) {
          int k = (int)(R::unif_rand() * (i + 1));
          std::swap(as_row[i], as_row[k]);
          std::swap(as_col[i], as_col[k]);
        }
        for (int idx = 0; idx < numActive; ++idx) {
          coord_update(as_row[idx], as_col[idx], J,
                       D, Psi, Rho, W, U, normD, diffD, S);
        }
        if (diffD <= normD * cdSweepTol) break;
      }
    }

    // compute logdetX on first pass
    if (fX >= 1e15) {
      mat chU;
      bool ok = chol(chU, Psi);
      if (!ok) chU = Psi; // fallback
      logdetX = 0.0;
      for (int i = 0; i < J; ++i) {
        double v = std::max(chU(i,i), 1e-6);
        logdetX += std::log(v);
      }
      logdetX *= 2.0;
      fX = (trSX + l1normX) - logdetX;
    }

    // trgradgD = tr((S - W) * D)
    double trgradgD = 0.0;
    for (int j = 0; j < J; ++j)
      trgradgD += (S(j,j) - W(j,j)) * D(j,j);
    for (int j = 1; j < J; ++j)
      for (int jj = 0; jj < j; ++jj)
        trgradgD += 2.0 * (S(j,jj) - W(j,jj)) * D(j,jj);

    double alpha = 1.0, l1normXD = 0.0, fX1prev = 1e15;

    for (int li = 0; li < max_lineiter; ++li) {
      double l1normX1 = 0.0, trSX1 = 0.0;
      mat Psi_next = Psi + alpha * D;
      compute_norms(J, Psi_next, S, Rho, l1normX1, trSX1);

      mat chPsi;
      bool ok = chol(chPsi, Psi_next);
      if (!ok) { alpha *= 0.5; continue; }

      double logdetX1 = 0.0;
      for (int l = 0; l < J; ++l) logdetX1 += std::log(chPsi(l,l));
      logdetX1 *= 2.0;
      fX1 = (trSX1 + l1normX1) - logdetX1;

      if (alpha == 1) l1normXD = l1normX1;

      if (fX1 <= fX + alpha * sigma * (trgradgD + l1normXD - l1normX) || normD == 0.0) {
        fXprev = fX; fX = fX1;
        l1normX = l1normX1; logdetX = logdetX1; trSX = trSX1;
        Psi = Psi + alpha * D;
        break;
      }
      if (fX1prev < fX1) {
        fXprev = fX;
        l1normX = l1normX1; logdetX = logdetX1; trSX = trSX1;
        Psi = Psi + alpha * D;
        break;
      }
      fX1prev = fX1;
      alpha *= 0.5;
    }

    W = inv(Psi);
    is_identity = false;

    bool converged = !(subgrad * alpha >= l1normX * eps &&
                       std::abs((fX - fXprev) / fX) >= EPS_INNER);
    if (converged) break;
    if (newton == max_iter_quic - 1)
      Rcpp::warning("hit max newton Iter within QUIC");
  }

  cube results(J, J, 2);
  results.slice(0) = Psi;
  results.slice(1) = W;
  return results;
}
