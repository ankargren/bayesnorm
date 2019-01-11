#ifndef BAYESNORM_MVN_CBM_H
#define BAYESNORM_MVN_CBM_H
inline arma::vec mvn_cbm(const arma::mat & Phi, const arma::vec & d, 
                         const arma::vec & alpha) {
  arma::uword n = Phi.n_rows;
  arma::uword p = Phi.n_cols;
  
  arma::mat U = Phi.t();
  U.each_col() %= d;
  arma::vec d_sqrt = sqrt(d);
  arma::mat I(n, n, arma::fill::eye);
  arma::vec u(p);
  u.imbue(norm_rand);
  arma::vec delta(n);
  delta.imbue(norm_rand);
  u %= d_sqrt;
  arma::vec v = Phi * u + delta;
  arma::vec w = arma::solve(Phi * U + I, (alpha - v));
  arma::vec theta = u + U * w;
  
  return theta;
}

#endif