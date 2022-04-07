#include <Rcpp.h>
using namespace Rcpp;


//count the unique elements in a vector
int countUnique2(NumericVector y);

// [[Rcpp::export]]
arma::vec poly_residuals(arma::vec yr, int m);  //use polynomial detrending and return residuals
arma::vec lm_c(arma::vec xs, arma::vec yr); // simple first order regression

arma::uvec seq_int(arma::uword length);