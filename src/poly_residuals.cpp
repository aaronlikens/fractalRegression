// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


arma::vec poly_residuals(arma::vec yr, int m){
  int rows = yr.n_rows;
  int cols = m + 1;
  arma::colvec coef(cols);
  
  arma::mat x(rows,cols);
  
  //convert data to an armadillo vector
  // arma::colvec y(yr.begin(), yr.size(), false);
  arma::colvec y = yr;
  
  //allocate memor for x and power of x vectors
  arma::colvec t1(rows);
  for ( int i = 0; i < rows; i++){
    t1(i) = i;
  }
  
  for ( int i = 0; i < cols; ++i){
    x.col(i) = arma::pow(t1,i);
  }
  
  coef = solve(x,y);
  
  //NumericVector coefs(coef.begin(),coef.end());
  
  arma::colvec resid = y-x*coef;
  
  //NumericVector Resid = NumericVector(resid.begin(),resid.end());
  
  return resid;
}
