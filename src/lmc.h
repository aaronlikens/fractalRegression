#include <Rcpp.h>
using namespace Rcpp;


//count the unique elements in a vector
int countUnique2(NumericVector y);

//use polynomial detrending and return residuals
//' Polynomial Residuals
//' Function that fits a polynomial and returns the residuals
//' @param yr is a real valued vector 
//' @param m is the detrending order
//' @export
// [[Rcpp::export]]
arma::vec poly_residuals(arma::vec yr, int m);  

// simple, fast first order regression
// Linear Model in c++
//' Simplef bivariate regression written in c++
//' @param xs a real valued column vector 
//' @param yr is a real valued column vector
//' @export
// [[Rcpp::export]]
arma::vec lm_c(arma::vec xs, arma::vec yr); 

// function to return a sequence of unsigned integers
//' Integer Sequence
//' Function that produces a sequence of integers from 1 to N
//' @param length is a positive integer that will produce a sequence from 1:length
//' @export
// [[Rcpp::export]]
arma::uvec seq_int(arma::uword length);

// function to return a sequence of unsigned integers defined on a range
//' Sequence of Integer ranges
//' Function that produces a sequece of integers that span a specific range
//' @param start is a positive integer and gives the smallest value in the sequence
//' @param stop is a positive integer and gives the largest value in a sequence
//' @export
// [[Rcpp::export]]
arma::uvec seq_int_range(arma::uword start, arma::uword stop);

// function to return detrended covariance
//' Detrended Covariance
//' Functional that returns the detrended covariance between two vectors
//' @param x a real valued column vector 
//' @param y is a real valued column vector
//' @param m is the detrending order
//' @export
// [[Rcpp::export]]
arma::vec detrend_cov(arma::vec x, arma::vec y, int m);