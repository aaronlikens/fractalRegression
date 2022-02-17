// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::export]]
arma::cx_vec make_complex(arma::vec x){
  arma::cx_vec cx_x(x,x);
  for (int i = 0; i < cx_x.n_rows; ++i){
    cx_x.row(i) = x[i];
  }
  
  return(cx_x);
  
}