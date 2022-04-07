// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


arma::uvec seq_int(arma::uword length){
  arma::uvec out(length);
  for (arma::uword i = 0; i < length; i++){
    out(i) = i;
  }
  return(out);
}