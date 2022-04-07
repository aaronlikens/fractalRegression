// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// function for doing the detrending. Returns the sum of squared residuals.
arma::vec detrend_cov(arma::vec x, arma::vec y, int m){
  int rows = x.n_elem;
  int cols = m + 1;
  arma::colvec coefx(cols);
  arma::colvec coefy(cols);
  
  arma::mat t(rows,cols);
  
  
  //allocate memor for x and power of x vectors
  arma::colvec t1(rows);
  for ( int i = 0; i < rows; i++){
    t1(i) = i+1;
  }
  //t1 = t1 - mean(t1);
  for ( int i = 0; i < cols; ++i){
    t.col(i) = arma::pow(t1,i);
  }
  
  // fit regression line
  coefx = solve(t,x);
  coefy = solve(t,y);
  
  // find residuals
  arma::colvec residx = x-t*coefx;
  arma::colvec residy = y-t*coefy;
  
  //square and sum the residuals
  double f2xy = 0;
  double f2x = 0;
  double f2y = 0;
  
  // find RMS and Covariance
  f2x = arma::accu(pow(residx,2))/(x.n_elem-1);
  f2y = arma::accu(pow(residy,2)/(x.n_elem-1));
  f2xy = arma::accu(residx % residy/(x.n_elem-1));
  
  arma::vec varCovar(3);
  varCovar[0] = f2x;
  varCovar[1] = f2y;
  varCovar[2] = f2xy;
  
  
  return varCovar;
}
