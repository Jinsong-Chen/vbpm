## Internal QUIC solver wrapper.
##
## `my_quic_cpp()` is the compiled RcppArmadillo solver (src/quic_impl.cpp,
## exported by Rcpp::compileAttributes into R/RcppExports.R). `.quic_rcpp()`
## keeps the interface used by vbfa()'s local-dependence branch: it returns a
## J x J x 2 array whose first slice is the residual precision Psi and whose
## second slice is its inverse W.
.quic_rcpp <- function(J, S, Rho, eps, max_iter_quic) {
  res <- my_quic_cpp(J, S, Rho, eps, max_iter_quic)
  out <- array(0, dim = c(J, J, 2))
  out[, , 1] <- res[, , 1]
  out[, , 2] <- res[, , 2]
  out
}
