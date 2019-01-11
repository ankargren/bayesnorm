#include <bayesnorm.h>

// [[Rcpp::export]]
arma::mat rmvn_cbm(const arma::mat & Phi, const arma::vec & d, 
                                       const arma::vec & alpha) {
    return mvn_cbm(Phi, d, alpha);
}
