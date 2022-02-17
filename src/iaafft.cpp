// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppArmadillo.h>
#include "fftw_fft.h"
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat iaafft_cpp(arma::vec signal, int N, int max_iter) {
    
  // define complex number i = sqrt(-1)
  std::complex<double> ii(0,1);
  arma::vec x = signal;
  int ln = x.n_elem;

  arma::cx_vec x_fft = fftw_fft(make_complex(x),'f');
  arma::vec amp = arma::abs(x_fft);
  arma::mat sgates(ln, N);
  sgates.zeros();
  
  for (uint n = 0; n < N; n++){
    arma::vec shuffle_x = arma::shuffle(x);
    sgates.col(n) = shuffle_x;
    }
  
  arma::uvec ind = sort_index(x);
  x = arma::sort(x);
//NOTE: UP to here, all processes match those in matlab
  for (uint n = 0; n < N; n++){
    arma::vec phase_x = arma::arg(fftw_fft(make_complex(sgates.col(n)),'f'));
    uint nn = 1;
    int converged = 0;
    arma::uvec ind_prev = ind;
    //check that routines has not covered or reached max iterations
    while (nn <= max_iter & converged == 0){
      arma::cx_vec shuffle_phase = amp%exp(phase_x*ii);
      arma::vec surrogate = arma::real(fftw_fft(shuffle_phase,'b'));
      arma::uvec indx = sort_index(surrogate);
      sgates.col(n) = arma::sort(surrogate);
      arma::uvec ind_new = sort_index(indx);
      sgates.col(n) = x.rows(ind_new);

      if (sum(ind_new == ind_prev)==ind_new.n_elem){
        converged = 1;
      }else{
        ind_prev = ind_new;
        nn +=1;
        if (nn == max_iter){

        }
        }
      phase_x = arma::arg(fftw_fft(make_complex(sgates.col(n)),'b'));
  
      }
    
    
    }
  

  
  return sgates;
}



