#ifndef BAYESNORM_MVN_H
#define BAYESNORM_MVN_H
inline arma::vec mvn_bcm(const arma::mat & Phi, const arma::vec & d, 
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

#ifndef BAYESNORM_MVN_RUE_H
#define BAYESNORM_MVN_RUE_H
inline arma::vec mvn_rue(const arma::mat & Phi, const arma::vec & d, 
                         const arma::vec & alpha) {
  
  arma::mat Q = Phi.t() * Phi;
  Q.diag() += pow(d, -1.0);
  arma::mat L = arma::chol(Q, "lower");
  arma::mat b = Phi.t() * alpha;
  arma::vec v = arma::solve(arma::trimatl(L), b);
  arma::vec mu = arma::solve(arma::trimatu(L.t()), v);
  arma::vec z(Phi.n_cols);
  z.imbue(norm_rand);
  arma::vec y = arma::solve(arma::trimatu(L.t()), z);
  arma::mat theta = mu + y;
  
  return theta;
}

#endif