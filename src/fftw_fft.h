#ifndef FFTW_FFT_H
#define FFTW_FFT_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
// #include "fftw3.h"


arma::cx_vec fftw_fft(arma::cx_colvec in_vec, char fft_type);
arma::cx_vec make_complex(arma::vec x);

#endif