## Internal wrapper for the compiled QUIC solver (src/quic_impl.cpp), used by
## vbfa()'s local-dependence branch. Returns a J x J x 2 array: slice 1 is the
## residual precision Psi, slice 2 is its inverse W.
.quic_rcpp <- function(J, S, Rho, eps, max_iter_quic)
  my_quic_cpp(J, S, Rho, eps, max_iter_quic)
