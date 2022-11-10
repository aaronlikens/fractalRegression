//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec win_sums(arma::vec x, unsigned int window);
arma::vec fitting(arma::vec x, arma::vec y);

// [[Rcpp::export]]
List mfdfa_cj(arma::vec Timeseries, arma::vec qValues, arma::uvec scales) {

  // initialise
  arma::uword nq = qValues.n_elem;
  arma::uword ns = scales.n_elem;
  arma::mat Ma = arma::zeros(nq,ns);
  arma::mat Mf = arma::zeros(nq,ns);
  arma::mat Md = arma::zeros(nq,ns);
  
  arma::vec muScale(scales.n_elem);
  for (unsigned int i = 0; i < scales.n_elem; i++){
    muScale(i) = -log10(pow(2, scales(i)));
  }
  double TimeseriesSummed = arma::accu(Timeseries);

  arma::vec alpha = arma::zeros(nq, 1);
  arma::vec falpha = arma::zeros(nq, 1);
  arma::vec Dq = arma::zeros(nq, 1);
  
  arma::vec Rsqr_alpha = arma::zeros(nq, 1);
  arma::vec Rsqr_falpha = arma::zeros(nq, 1);
  arma::vec Rsqr_Dq = arma::zeros(nq, 1);
  for (unsigned int i = 0; i < nq; i++){
    
    double q = qValues(i);
    
    for (unsigned int j = 0; j < ns; ++j){

      // determine how many windows we will have at this scale
      unsigned int window = pow(2, scales(j));

      //break the time series into windows & sum
      arma::vec ps = win_sums(Timeseries, window);
      arma::vec p = ps/TimeseriesSummed;
  
      double Nor = arma::accu(arma::pow(p, q));

      //calculation of Md
      Md(i, j) = log10(Nor); //%not accounting for q between 0 and 1
      if (q <=1 & q > 0){
        Md(i, j) = arma::accu(p % arma::log10(p))/Nor;
      }
      
      // Ma & Mf
      arma::vec mu = arma::pow(p, q)/Nor;
      Ma(i,j) = arma::accu(mu % log10(p));
      Mf(i,j) = arma::accu(mu % log10(mu));

    }
    
    // regression part
    arma::vec b_Ma = fitting(muScale, Ma.row(i).t());
    arma::vec b_Mf = fitting(muScale, Mf.row(i).t());
    arma::vec b_Md = fitting(muScale, Md.row(i).t());
    
    alpha(i) = b_Ma(0);
    falpha(i) = b_Mf(0);
    Dq(i) = b_Md(0)/(q - 1);
    if (q <= 1 & q > 0){
      Dq(i) = b_Md(0);
    }

 
    Rsqr_alpha(i) = b_Ma(1);
    Rsqr_falpha(i) = b_Mf(1);
    Rsqr_Dq(i) = b_Md(1); 

  }
  
  
  return List::create(Named("x") = Timeseries,
                      Named("alpha") = alpha,
                      Named("falpha") = falpha,
                      Named("Dq") = Dq,
                      Named("Rsqr_alpha") = Rsqr_alpha,
                      Named("Rsqsr_falpha") = Rsqr_falpha,
                      Named("Rsqr_Dq") = Rsqr_Dq,
                      Named("muScale") = muScale,
                      Named("Md") = Md,
                      Named("Ma") = Ma,
                      Named("Md") = Mf,
                      Named("q") = qValues,
                      Named("scales") = scales);
}


// divide time series into windows and compute the sum of each window;
// [[Rcpp::export]]
arma::vec win_sums(arma::vec x, unsigned int window){
  arma::uword winsize = floor(x.n_elem/window);
  arma::uword nwindows = window;
  arma::vec out(nwindows);
  arma::uword start= 0;
  arma::uword stop = winsize-1;
  
  for (unsigned int i = 0; i < nwindows; i++){
    out(i) = arma::accu(x.rows(start,stop));
    start += winsize;
    stop += winsize;
  }
  return out;
}

// [[Rcpp::export]]
arma::vec fitting(arma::vec x, arma::vec y){
  
  // CALCULATING THE COEFFICIENTS
  arma:: mat X = arma::ones(x.n_elem,2);
  X.col(1) = x;
  arma::vec B = arma::solve(X, y);
  arma::vec yF = X*B;
  
  
  // CALCULATING THE R2
  double R2 = 1 - arma::accu(arma::pow(y - yF, 2))/arma::accu(arma::pow(y - mean(y),2));
  
  B(0) = B(1);
  B(1) = R2;
  return B;
  
  
}