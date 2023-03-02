// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// function for lagging variables
arma::mat lagn(arma::vec x, int lags);

// functions for detrending and computing scale-wise covariance
arma::vec polyfit(arma::vec x, int order);
arma::mat detrend_var(arma::mat x, int order);

//' Multiscale Lagged Regression Analysis
 //'
 //' Fast function for computing multiscale lagged regression analysis (MLRA) on long time series. Combining DFA with ordinary least square regression, MLRA
 //' is a form of fractal regression that can be used to estimate asymmetric and multiscale regression coefficients between two variables at different time-scales and temporal lags. 
 //'
 //' @param x A real valued vector (i.e., time series data) to be analyzed.
 //' @param y A real valued vector (i.e., time series data) to be analyzed.
 //' @param order is an integer indicating the polynomial order used for 
 //' detrending the local windows (e.g, 1 = linear, 2 = quadratic, etc.). There 
 //' is a not pre-determined limit on the order of the polynomial order but the 
 //' user should avoid using a large polynomial on small windows. This can result
 //' in overfitting and non-meaningful estimates. 
 //' @param scales An integer vector of scales over which to compute correlation. 
 //' Unlike univariate DFA, MRA does not require that scales be in log units.
 //' Scale intervals can be sequential, for example, when the analysis is 
 //' exploratory and no a priori hypotheses have been made about the scale of 
 //' correlation. A small subset of targeted scales may also be investigated 
 //' where scale-specific research questions exist. We have found that windows
 //' smaller than say 8 observations create stability problems due to 
 //' overfitting. This is espcially when the order of the fitting polynomial is 
 //' large.
 //' @param lags An integer indicating the maximum number of lags to include in the analysis.
 //' @param direction A character string indicating a positive ('p') or negative ('n') lag.
 //' @import Rcpp
 //' @useDynLib fractalRegression
 //' @export
 //'
 //' @details Mathematical treatment of the MLRA algorithm and its performance can be found in Kristoufek (2015) and Likens et al. (2019).
 //'
 //' Use of the direction parameter specifies whether the scale-wise \eqn{\beta} coefficients for positive or negative lags will be estimated.  
 //'
 //' Note that under conditions with linear and quadratic trends, Likens et al. (2019) found that there was a systematic positive bias in the \eqn{\beta} estimates for larger scales.
 //' Using a polynomial detrending order of 2 or greater was shown to attenuate this bias. 
 //'
 //' @return The object returned from the mlra() function is a list containing \code{betas} the \eqn{\beta} coefficients for each lag at each of the scales. 
 //'
 //'
 //' @references
 //'
 //' Kristoufek, L. (2015). Detrended fluctuation analysis as a regression framework: Estimating dependence at different scales. Physical Review E, 91(2), 022802.
 //'
 //' Likens, A. D., Amazeen, P. G., West, S. G., & Gibbons, C. T. (2019). Statistical properties of Multiscale Regression Analysis: Simulation and application to human postural control. Physica A: Statistical Mechanics and its Applications, 532, 121580.
 //'
 //' @examples
 //'
 //' 
 //' # Here is a simple example for running MLRA using a white noise and pink noise time series.
 //' # For more detailed examples, see the vignette. 
 //' 
 //' noise <- rnorm(5000)
 //' 
 //' pink.noise <- fgn_sim(n = 5000, H = 0.9)
 //'
 //' scales <- logscale(scale_min = 10, scale_max = 1250, scale_ratio = 1.1)
 //' 
 //' mlra.out <- mlra(
 //'     x = noise, 
 //'     y = pink.noise, 
 //'     order = 1, 
 //'     scales = scales, 
 //'     lags = 100, direction = 'p')
 //' 
 //'
 //'
 //[[Rcpp::export]]
 List mlra(arma::vec x, arma::vec y, int order, IntegerVector scales,
            int lags, char direction){
   
   // should be mat(N-l-1, l + 1), where l = number of lags
   arma::mat lx(x.n_elem - lags, lags + 1, arma::fill::zeros);
   
   // should be mat(N-l-1, l + 1), where l = number of lags
   arma::mat ly(x.n_elem - lags, lags + 1, arma::fill::zeros); 
   
   // should be mat(N-l-1, l + 2), where l = number of lags
   arma::mat lyx(x.n_elem - lags, lags + 2, arma::fill::zeros); 
   
   unsigned int n_scales = scales.size();
   NumericVector resid(n_scales);
   arma::field<arma::mat> rmse(n_scales, 1);
   arma::mat betas(lags + 1, n_scales);
   
   arma::vec X = cumsum(x-mean(x));
   arma::vec Y = cumsum(y-mean(y));
   
   
   // negative
   if (direction == 'n'){
     lx = lagn(X, lags);
     ly = lagn(Y, lags);
     lyx.resize(ly.n_rows, ly.n_cols + 1);
     lyx.col(lyx.n_cols-1) = ly.col(0);
     for (unsigned int i = 0; i < lyx.n_cols-1; ++i){
       lyx.col(i) = lx.col(i);
     }
     // positive
   }else if (direction == 'p'){
     ly = lagn(Y, lags);
     lx = lagn(X, lags);
     lyx.resize(ly.n_rows, ly.n_cols + 1);
     lyx.col(lyx.n_cols-1) = lx.col(0);
     // lyx.cols(1,lyx.n_cols) = ly;
     for (unsigned int i = 0; i < lyx.n_cols-1; ++i){
       lyx.col(i) = ly.col(i);
     }
     
   }
   unsigned int p = lyx.n_cols;
   unsigned int T = lyx.n_rows;
   // allocate memory to store covariance matrices
   arma::vec sxy(p-1, arma::fill::zeros);
   arma::mat sxx(p-1, p-1, arma::fill::zeros);
   arma::mat s(p, p, arma::fill::zeros);
   for (unsigned int i = 0; i < n_scales; ++i){
     unsigned int start = 0; 
     unsigned int stop = scales[i];
     arma::mat varcov(p,p, arma::fill::zeros);
     int count = 0; 
     while (stop < T-1){
       
       varcov += detrend_var(lyx(arma::span(start,stop), arma::span::all), 
                             order);
       start = start + scales[i];
       stop = stop + scales[i];
       count += 1;
       //Rcout << "Start = " << start <<". Stop = " << stop << "\n";
     }
     rmse(i,0) = varcov/(count*scales[i]);
     s = rmse(i,0);
     
     sxy = s(arma::span(0,s.n_rows-2), s.n_cols-1);
     sxx = s(arma::span(0,s.n_rows-2), arma::span(0,s.n_cols-2));
     Rcout << sxy.n_rows << " , " << sxx.n_cols << "\n";
     // betas.col(i) = solve(sxx, sxy);
     // calculate scale-wise regression coefficients
     betas.col(i) = arma::inv(sxx)*sxy;
   }
   
   return List::create(Named("betas") = betas);
 }
 
 //[[Rcpp::export]]
 arma::mat lagn(arma::vec x, int lags) {
   
   arma::mat lx(x.n_elem - lags, lags+1, arma::fill::zeros);
   int L = lags;
   unsigned int start = L;
   unsigned int stop = x.n_elem-1;
   
   // negative lag
   for (int i = 0; i <= L; ++i){
     // lx.col(i) = x(arma::span(L-i-1, x.n_elem-i-1));
     lx.col(i) = x(arma::span(start, stop));
     start -= 1;
     stop -= 1;
     
   }
   
   return lx;
 }
 
 //[[Rcpp::export]]
 arma::vec polyfit(arma::vec x, int order){
   unsigned int cols = order + 1;
   unsigned int rows = x.n_elem;
   
   arma::colvec coefx(cols);
   arma::mat t(rows,cols);
   
   // create time vectors up to specified polynomial order
   arma::colvec t1(rows);
   t1 = arma::regspace(1, 1, rows);
   for ( unsigned int i = 0; i < cols; ++i){
     t.col(i) = arma::pow(t1,i);
   }
   
   coefx = solve(t,x);
   arma::vec residx = x-t*coefx;
   return residx;
 }
 
 //[[Rcpp::export]]
 arma::mat detrend_var(arma::mat X, int order){
   
   arma::mat d(X.n_rows, X.n_cols, arma::fill::zeros);
   
   for (unsigned int i = 0; i < X.n_cols; ++i){
     d.col(i) = polyfit(X.col(i), order);
   }
   
   arma::mat varcov = cov(d);
   // arma::mat varcov = inv(d.t()*d);
   return varcov;
   
 }
