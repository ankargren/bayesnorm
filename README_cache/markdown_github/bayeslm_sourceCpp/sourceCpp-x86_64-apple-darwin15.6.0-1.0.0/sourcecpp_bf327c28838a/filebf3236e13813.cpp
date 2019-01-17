#include <RcppArmadillo.h>
#include <mvnorm.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::export]]
arma::mat bayeslm(const arma::vec& y, const arma::mat& x,
                   const int iters = 1000) {
  int p = x.n_cols;
  arma::vec d = arma::vec(p, arma::fill::ones); // prior variance is 1
  arma::mat xtx, Sigma, mu;
  
  // Storage
  arma::mat beta_draws(iters, p); // Object to store beta draws in
  for ( int iter = 0; iter < iters; ++iter ) {
    xtx = x.t() * x; // X'X
    xtx.diag() += arma::pow(d, -1.0); // add D^{-1}
    
    Sigma = xtx.i(); // the inverse is Sigma
    mu = Sigma * x.t() * y; // compute mu
    
    beta_draws.row(iter) = rmvnorm(1, mu, Sigma);
  }
  return beta_draws;
}


#include <Rcpp.h>
// bayeslm
arma::mat bayeslm(const arma::vec& y, const arma::mat& x, const int iters);
RcppExport SEXP sourceCpp_3_bayeslm(SEXP ySEXP, SEXP xSEXP, SEXP itersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const int >::type iters(itersSEXP);
    rcpp_result_gen = Rcpp::wrap(bayeslm(y, x, iters));
    return rcpp_result_gen;
END_RCPP
}
