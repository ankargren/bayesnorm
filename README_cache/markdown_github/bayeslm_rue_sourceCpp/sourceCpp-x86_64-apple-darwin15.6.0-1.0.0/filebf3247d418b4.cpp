#include <RcppArmadillo.h>
#include <bayesnorm.h>
// [[Rcpp::depends(RcppArmadillo, bayesnorm)]]
// [[Rcpp::export]]
arma::mat bayeslm_rue(const arma::vec& y, const arma::mat& x,
                       const int iters = 1000) {
  int p = x.n_cols;
  arma::vec d = arma::vec(p, arma::fill::ones);
  arma::mat beta_draws(p, iters);
  for ( int iter = 0; iter < iters; ++iter ) {
    beta_draws.col(iter) = mvn_rue(x, d, y);
  }
  return beta_draws;
}
