//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "lmc.h"
using namespace Rcpp;

//' Multifractal Detrended Fluctuation Analysis 2
//'
//' Fast function for computing multifractal detrended fluctuation analysis 
//' (MF-DFA), a widely used method for estimating the family of long-range 
//' temporal correlations or scaling exponents in time series data. 
//' MF-DFA is also a form of multifractal analysis that indicates the degree 
//' of interaction across temporal scales.
//' 
//' @param x A real valued vector (i.e., time series data) to be analyzed. 
//' @param order is an integer indicating the polynomial order used for 
//' detrending the local windows (e.g, 1 = linear, 2 = quadratic, etc.). There 
//' is not pre-determined limit on the order of the polynomial order but the 
//' user should avoid using a large polynomial on small windows. This can result
//' in overfitting and non-meaningful estimates. 
//' @param scales An integer valued vector indicating the scales one wishes to resolve
//' in the analysis. Best practice is to use scales which are evenly spaced in 
//' the logarithmic domain e.g., \code{scales = 2^(4:(N/4))}, where N is the length of the
//' time series. Other, logarithmic bases may also be used to give finer 
//' resolution of scales while maintaining ~= spacing in the log domain e.g, 
//' \code{scales = unique(floor(1.1^(30:(N/4))))}. Note that fractional bases may 
//' produce duplicate values after the necessary floor function.
//' @import Rcpp
//' @useDynLib fractalRegression
//' @export
//' 
//' @details Details of the algorithm are specified in detail in Kantelhardt et al. (2001; 2002) and visualized nicely in Kelty-Stephen et al. (2016).
//' 
//' Selecting the range of values for q is important. Note that MF-DFA estimates for q = 2 are equivalent to DFA. Larger values of q (q > 2) emphasize larger residuals and smaller values of q
//' (q < 2) emphasis smaller residuals (Kelty-Stephen et al., 2016). For most biomedical signals such as physiological and kinematic, a q range of -5 to 5 is common (Ihlen, 2010). However, in some cases, 
//' such as when time series are short (< 3000), it can be appropriate to limit the range of q to positive only. Kelty-Stephen et al. (2016) recommend a 
//' positive q range of 0.5 to 10 with an increment of 0.5. 
//'
//' While it is common to use only linear detrending with DFA and MF-DFA, it is important to inspect the trends in the data to determine
//' if it would be more appropriate to use a higher order polynomial for detrending, and/or compare the DFA and MF-DFA output for different polynomial orders (see Ihlen, 2012; Kantelhardt et al., 2001).
//' 
//' General recommendations for choosing the min and max scale are a scale_min = 10 and scale_max = (N/4), where N is the number of observations.
//' See Eke et al. (2002), Gulich and Zunino (2014), Ihlen (2012), and  for additional considerations and information on choosing the correct parameters. 
//'
//' @return The output of the algorithm is a list that includes:
//' \itemize{ 
//'  \item \code{log_scale} The log scales used for the analysis
//'  \item \code{log_fq} The log of the fluctuation functions for each scale and q 
//'  \item \code{Hq} The q-order Hurst exponent (generalized Hurst exponent)
//'  \item \code{Tau} The q-order mass exponent
//'  \item \code{q} The q-order statistical moments
//'  \item \code{h} The q-order singularity exponent
//'  \item \code{Dh} The dimension of the q-order singularity exponent
//'}
//'
//' @references 
//'
//' Ihlen, E. A. F. (2012). Introduction to Multifractal Detrended Fluctuation Analysis in Matlab. Frontiers in Physiology, 3. https://doi.org/10.3389/fphys.2012.00141
//'
//' Kantelhardt, J. W., Koscielny-Bunde, E., Rego, H. H., Havlin, S., & Bunde, A. (2001). Detecting long-range correlations with detrended fluctuation analysis. Physica A: Statistical Mechanics and its Applications, 295(3-4), 441-454.
//' 
//' Kantelhardt, J. W., Zschiegner, S. A., Koscielny-Bunde, E., Havlin, S., Bunde, A., & Stanley, H. E. (2002). Multifractal detrended fluctuation analysis of nonstationary time series. Physica A: Statistical Mechanics and its Applications, 316(1-4), 87-114.
//'
//' Kelty-Stephen, D. G., Palatinus, K., Saltzman, E., & Dixon, J. A. (2013). A Tutorial on Multifractality, Cascades, and Interactivity for Empirical Time Series in Ecological Science. Ecological Psychology, 25(1), 1-62. https://doi.org/10.1080/10407413.2013.753804
//'
//' Kelty-Stephen, D. G., Stirling, L. A., & Lipsitz, L. A. (2016). Multifractal temporal correlations in circle-tracing behaviors are associated with the executive function of rule-switching assessed by the Trail Making Test. Psychological Assessment, 28(2), 171-180. https://doi.org/10.1037/pas0000177
//'
//' @examples
//'
//' 
//' 
//' noise <- rnorm(5000)
//' 
//' scales <- c(16,32,64,128,256,512,1024)
//'
//' mf.dfa.white.out <- mfdfa(
//'     x = noise, q = c(-5:5), 
//'     order = 1, 
//'     scales = scales, 
//'     scale_ratio = 2) 
//'  
//' pink.noise <- fgn_sim(n = 5000, H = 0.9)
//' 
//' mf.dfa.pink.out <- mfdfa(
//'     x = pink.noise, 
//'     q = c(-5:5), 
//'     order = 1, 
//'     scales = scales, 
//'     scale_ratio = 2)
//'
//' 
//' 
// [[Rcpp::export]]
List mfdirect(arma::vec x, int order, arma::uvec scales) {

  arma::vec Fq0(scales.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < scales.n_elem; i++){
    arma::uword nbins = floor(x.n_elem/scales(i) );
    arma::vec rms0(nbins, arma::fill::zeros);
    arma::uvec indx = seq_int_range(0, scales(i) - 1);
    
    for (arma::uword j = 0; j <  nbins; j++ ){
      arma::vec temp = x.rows(indx);
      arma::vec res = poly_residuals(temp, order);
      res = pow(res,2);
      rms0(j) = sqrt(mean(res));
      indx = indx + scales(i);
    }
    rms0 = pow(rms0,2);
    Fq0(i) = exp(0.5*mean(log(rms0)));
  }
  
  
  arma::uword halfmax = floor((double)max(scales)/2);
  arma::uvec time_index = seq_int_range(halfmax - 1, x.n_elem - halfmax - 1);


  arma::vec F(scales.n_elem);
  arma::field<arma::vec> RMS(scales.n_elem, 1);
  for (arma::uword i = 0; i < scales.n_elem; i++){
    
    arma::uword halfseg = floor((double)scales(i)/2);
    arma::uword nbins = x.n_elem - 2*halfmax - 1;
    arma::vec rms(nbins + 1);
    arma::uword counter = 0;
    for (arma::uword j = halfmax - 1; j < x.n_elem - halfmax - 1; j ++){  
      // Rcout << "Outer Loop " << i << " out of " << scales.n_elem << "\n";
      // Rcout << "halfseg is " << halfseg << "\n";
      // Rcout << "Inner Loop " << j << " out of " << nbins << "\n";
      arma::uvec indx = seq_int_range(j - halfseg - 1, j + halfseg - 1);
      // Rcout << "First and Last indices are " << indx(0) << " and " <<
      //   indx(indx.n_elem-1) << ". x has " << x.n_rows << "rows.\n";
      arma::vec temp = x.rows(indx);
      // Rcout << "I indexed the rows.\n";
      arma::vec res = poly_residuals(temp, order);
      // Rcout << "I did the regression.\n";
      res = pow(res,2);
      // Rcout << "j is " << j << " and rms has " << rms.n_elem << "rows.\n";
      rms(counter) = sqrt(mean(res));
      counter += 1;
      // Rcout << "I wrote to the rms variable.\n";
      
    }
    RMS(i) = rms;
    rms = pow(rms,2);
    F(i) = exp(0.5*mean(log(rms)));
  }
  // Rcout << "Madw it thorugh the second loop \n";
  arma::vec logs = arma::log2(arma::conv_to<arma::vec>::from(scales));
  arma::vec logfq0 = log2(Fq0);
  arma::vec C = lm_c(logs, logfq0);
  
  // create augmented matrix
  arma::mat aug_logs(logs.n_elem, 2, arma::fill::ones);
  aug_logs.col(1) = logs;
  
  // get predicted values
  arma::vec regfit =  aug_logs*C;
  
  double Hq0=C(1);
  arma::uword maxL = time_index(time_index.n_elem-1);
  arma::mat Ht(time_index.n_elem - 1, scales.n_elem);
  // Rcout << "about to enter the loop\n";
  for (arma::uword i = 0; i < scales.n_elem; ++i){
    // Rcout << "Loop " << i << " out of " << scales.n_elem << "\n";
    // Rcout << "First and Last indices are " << time_index(0) << " and " <<
    //   time_index(time_index.n_elem-1) << ".\n";
    // Rcout << "RMS has " << RMS(i).n_rows << " rows.\n";
    arma::vec rmst = RMS(i);
    // arma::vec rmst = RMS(i).rows(time_index);
    Rcout << "got the RMS\n";
    arma::vec resrms = regfit(i) - log2(rmst);
    double logscale = log2(maxL) - log2(scales(i));
    // Rcout << "Is Ht the problem?" << "\n";
    Ht.col(i) = (resrms/logscale) + Hq0;
    // Rcout << "maybe not\n";
  }
  // Rcout << "Made it thorugh the third loop \n";
  //NOTE: should this be a floor?
  Rcout << "1Here.\n";
  double bin_numb = round(sqrt(Ht.n_elem));
  Rcout << "2Here.\n";
  // Ht.resize(Ht.n_elem);
  Rcout << "3Here.\n";
  double range_Ht = Ht.max() - Ht.min();
  Rcout << "4Here.\n";
  double binsize = Ht.n_elem/range_Ht;
  Rcout << "5Here.\n";
  arma::vec Htbin(bin_numb, arma::fill::zeros);
  Rcout << "6Here.\n";
  Htbin(0) = Ht.min();
  Rcout << "7Here.\n";
  for (arma::uword i = 1; i < Htbin.n_elem; ++i){
    Htbin(i) = Htbin(i - 1) + binsize;
  }
  Rcout << "8Here.\n";
  arma::vec Htrow = Ht;
  Htrow.resize(Htrow.n_elem);
  arma::uvec freq = hist(Htrow, Htbin);
  
  
  Rcout << "9Here.\n";
  arma::vec Ph = arma::conv_to<arma::vec>::from(freq)/sum(freq);
  Rcout << "10Here.\n";
  arma::vec Ph_norm = Ph/max(Ph);
  Rcout << "11Here.\n";
  arma::vec Dh = 1- (log(Ph_norm))/log(mean(diff(Htbin)));
  Rcout << "12Here.\n";
  
  return List::create(Named("x") = x,
                      Named("scales") = scales,
                      Named("Ht") = Ht,
                      Named("htbin") = Htbin,
                      Named("Ph") = Ph,
                      Named("Ph_norm") = Ph_norm,
                      Named("Dh") = Dh);
}


