#include <Rcpp.h>
using namespace Rcpp;


//count the unique elements in a vector
int countUnique2(NumericVector y);

//use polynomial detrending and return residuals
// [[Rcpp::export]]
arma::vec poly_residuals(arma::vec yr, int m);  

// simple, fast first order regression
// [[Rcpp::export]]
arma::vec lm_c(arma::vec xs, arma::vec yr); 

// function to return a sequence of unsigned integers
// [[Rcpp::export]]
arma::uvec seq_int(arma::uword length);

// function to return detrended covariance
// [[Rcpp::export]]
arma::vec detrend_cov(arma::vec x, arma::vec y, int m);