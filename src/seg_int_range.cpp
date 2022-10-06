// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


arma::uvec seq_int_range(arma::uword start, arma::uword stop){
  arma::uvec out(stop - start + 1);
  arma::uword counter = 0;
  for (arma::uword i = start; i <= stop; i++){
    out(counter) = i;
    counter += 1;
  }
  return(out);
}