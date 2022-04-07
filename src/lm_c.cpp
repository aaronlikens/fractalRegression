// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//compute simple linear slope and intercept
arma::vec lm_c(arma::vec xs, arma::vec yr){
  arma::vec coefs;
  int n = yr.size();
  
  // create augmented matrix
  arma::mat augmentedX(n, 2, arma::fill::ones);
  augmentedX.col(1) = xs;
  
  //compute regression coefficients
  arma::colvec coef = arma::solve(augmentedX, yr, arma::solve_opts::fast);
  
  return coef;
}
