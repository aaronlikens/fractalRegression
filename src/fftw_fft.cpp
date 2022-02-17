// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
#include "fftw3.h"


//' Fastest Fourier Transform in the west
//' A simple 1D implementation of the FFTW and its inverse
//' @param in_vec a real or complex valued vector.
//' @param type is a character with two values, 'f' = forward, 'i' = inverse
//' @import Rcpp
//' @useDynLib fractalRegression
//' @export
//' @return returns a complex valued 1D FFT
//[[Rcpp::export]]
arma::cx_vec fftw_fft(arma::cx_colvec in_vec, char fft_type){
  int matrix_size = (int)in_vec.n_rows;
  arma::cx_colvec out_vec(in_vec.n_elem, arma::fill::zeros);
  
  fftw_plan p;
  if (fft_type == 'f'){
    p = fftw_plan_dft_1d(matrix_size,
                         // invec,
                         (fftw_complex*) in_vec.memptr(),
                         // outvec,
                         reinterpret_cast<fftw_complex*>(out_vec.memptr()),
                         FFTW_FORWARD,
                         FFTW_ESTIMATE);
    
  }else{
    p = fftw_plan_dft_1d(matrix_size,
                         // invec,
                         (fftw_complex*) in_vec.memptr(),
                         // outvec,
                         reinterpret_cast<fftw_complex*>(out_vec.memptr()),
                         FFTW_BACKWARD,
                         FFTW_ESTIMATE);
    
  }
  
  fftw_execute(p);
  fftw_destroy_plan(p);
  if (fft_type == 'f'){
    return out_vec;  
  }else {
    return out_vec/matrix_size;
  }
  
}


