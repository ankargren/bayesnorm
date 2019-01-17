#include <RcppArmadillo.h>
#include <bayesnorm.h>
// [[Rcpp::depends(RcppArmadillo, bayesnorm)]]
// [[Rcpp::export]]
arma::mat bayeslm_bcm(const arma::vec& y, const arma::mat& x,
                       const int iters = 1000) {
  int p = x.n_cols;
  arma::vec d = arma::vec(p, arma::fill::ones);
  arma::mat beta_draws(p, iters);
  for ( int iter = 0; iter < iters; ++iter ) {
    beta_draws.col(iter) = mvn_bcm(x, d, y);
  }
  return beta_draws;
}


#include <Rcpp.h>
// bayeslm_bcm
arma::mat bayeslm_bcm(const arma::vec& y, const arma::mat& x, const int iters);
RcppExport SEXP sourceCpp_3_bayeslm_bcm(SEXP ySEXP, SEXP xSEXP, SEXP itersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type iters(itersSEXP);
    rcpp_result_gen = Rcpp::wrap(bayeslm_bcm(y, x, iters));
    return rcpp_result_gen;
END_RCPP
}
