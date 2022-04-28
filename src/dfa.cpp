// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "lmc.h"
using namespace Rcpp;


//' Detrended Fluctuation Analysis
//' 
//' Fast function for computing detrended fluctuation analysis (DFA), a widely used method for estimating long-range temporal correlations in time series data. 
//' DFA is also a form of mono-fractal analysis that indicates the degree of self-similarity across temporal scales.
//' 
//' @param x A real valued vector (i.e., time series data) to be analyzed. 
//' @param order An integer indicating the polynomial order used for 
//' detrending the local windows (e.g, 1 = linear, 2 = quadratic, etc.). There 
//' is not a pre-determined limit on the order of the polynomial order but the 
//' user should avoid using a large polynomial on small windows. This can result
//' in overfitting and non-meaningful estimates. 
//' @param verbose If the value of verbose = 1, then a list object is returned 
//' that includes: \code{log_scales} the log of all included scales, 
//' \code{log_rms} the log root mean square error (RMS) per scale, 
//' and \code{alpha} the overall \eqn{\alpha} estimate. If the value of 
//' verbose = 0, then a list containing only `alpha` will be returned.
//' @param scales An integer valued vector indicating the scales one wishes to resolve
//' in the analysis. Best practice is to use scales which are evenly spaced in 
//' the logarithmic domain e.g., \code{scales = 2^(4:(N/4))}, where N is the length of the
//' time series. Other, logarithmic bases may also be used to give finer 
//' resolution of scales while maintaining ~= spacing in the log domain e.g, 
//' \code{scales = unique(floor(1.1^(30:(N/4))))}. Note that fractional bases may 
//' produce duplicate values after the necessary floor function.
//' @param scale_ratio A scaling factor by which successive window sizes were 
//' were created. The default is 2 but should be addressed according to how 
//' scales were generated for example using \code{logscale(16, 100, 1.1)}, 
//' where 1.1 is the scale ratio.
//' @import Rcpp
//' @useDynLib fractalRegression
//' @export
//' 
//' @details Details of the algorithm are specified in detail in Peng et al. (1994) and visualized nicely in Kelty-Stephen et al. (2016).
//' The output of the algorithm is an \eqn{\alpha} (alpha) estimate which is a generalization of the Hurst Exponent. Conventional interpretation of \eqn{\alpha} is:
//' \itemize{
//'  \item \eqn{\alpha < 0.5 =} anti-correlated
//'  \item \eqn{\alpha ~= 0.5 =} uncorrelated, white noise
//'  \item \eqn{\alpha > 0.5 =} temporally correlated
//'  \item \eqn{\alpha ~= 1 =} 1/f-noise, pink noise
//'  \item \eqn{\alpha > 1 =} non-stationary and unbounded
//'  \item \eqn{\alpha ~= 1.5 =} fractional brownian motion
//' } 
//' 
//' We recommend a few points of consideration here in using this function. One is to be sure to 
//' verify there are not cross-over points in the logScale-logFluctuation plots (Peng et al., 1995; Perakakis et al ., 2009). Cross-over points 
//' (or a visible change in the slope as a function of of scale) indicate that a mono-fractal characterization 
//' does not sufficiently characterize the data. If cross-over points are evident, we recommend proceeding to using the mfdfa() to estimate the multi-fractal
//' fluctuation dynamics across scales.
//' 
//' While it is common to use only linear detrending with DFA, it is important to inspect the trends in the data to determine
//' if it would be more appropriate to use a higher order polynomial for detrending, and/or compare the DFA output for different polynomial orders (see Kantelhardt et al., 2001).
//' 
//' General recommendations for choosing the min and max scale are an sc_min = 10 and sc_max = (N/4), where N is the number of observations.
//' See Eke et al. (2002) and Gulich and Zunino (2014) for additional considerations. 
//' 
//' @return The object returned can take the following forms:
//' \itemize{ 
//'  \item If the value of verbose = 1, then a list object is returned that includes: \code{log_scales}
//' the log of all included scales, \code{log_rms} the log root mean square error (RMS) per scale, and \code{alpha} the overall \eqn{\alpha} estimate.
//'  \item If the value of verbose = 0, then a list containing only `alpha` the estimated scaling exponent \eqn{\alpha} will be returned.
//' }
//' @references 
//' 
//' Eke, A., Herman, P., Kocsis, L., & Kozak, L. R. (2002). Fractal characterization of complexity in temporal physiological signals. Physiological measurement, 23(1), R1-R38.
//' 
//' Gulich, D., & Zunino, L. (2014). A criterion for the determination of optimal scaling ranges in DFA and MF-DFA. Physica A: Statistical Mechanics and its Applications, 397, 17-30.
//'
//' Kantelhardt, J. W., Koscielny-Bunde, E., Rego, H. H., Havlin, S., & Bunde, A. (2001). Detecting long-range correlations with detrended fluctuation analysis. Physica A: Statistical Mechanics and its Applications, 295(3-4), 441-454.
//' 
//' Kelty-Stephen, D. G., Stirling, L. A., & Lipsitz, L. A. (2016). Multifractal temporal correlations in circle-tracing behaviors are associated with the executive function of rule-switching assessed by the Trail Making Test. Psychological assessment, 28(2), 171-180.
//' 
//' Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, and Goldberger AL (1994), Mosaic organization of DNA nucleotides, Physical Review E, 49, 1685-1689.
//' 
//' Peng C-K, Havlin S, Stanley HE, and Goldberger AL (1995), Quantification of scaling exponents and crossover phenomena in nonstationary heartbeat time series, Chaos, 5, 82-87.
//' 
//' Perakakis, P., Taylor, M., Martinez-Nieto, E., Revithi, I., & Vila, J. (2009). Breathing frequency bias in fractal analysis of heart rate variability. Biological psychology, 82(1), 82-88.
//' 
//' @examples
//' 
//' 
//' 
//' noise <- rnorm(5000)
//'
//' scales <- c(16,32,64,128,256,512,1024)
//
//' dfa.noise.out <- dfa(
//'     x = noise, 
//'     order = 1, 
//'     verbose = 1, 
//'     scales = scales,
//'     scale_ratio = 2)
//' 
//' pink.noise <- fgn_sim(n = 5000, H = 0.9)
//' 
//' dfa.pink.out <- dfa(
//'     x = pink.noise, 
//'     order = 1, 
//'     verbose = 1, 
//'     scales = scales, 
//'     scale_ratio = 2)
//' 
//' anticorr.noise <- fgn_sim(n = 5000, H = 0.25)
//' 
//' dfa.anticorr.out <- dfa(
//'     x = anticorr.noise, 
//'     order = 1, 
//'     verbose = 1, 
//'     scales = scales, 
//'     scale_ratio = 2)
//'   
//' 
//' 
//' 
// [[Rcpp::export]]
List dfa(arma::vec x, int order, arma::uword verbose, 
             arma::uvec scales, double scale_ratio = 2){
    
    double len = x.n_elem;
    arma::uword number_of_scales = scales.n_elem;
    arma::vec resid(number_of_scales);
    arma::vec X = cumsum(x-mean(x));
    
    // do the detrending and return the RMSE for each of the ith scales
    arma::vec RMS(number_of_scales);
    
    for (arma::uword i = 0; i < number_of_scales; ++i){
        arma::uword window = scales[i];
        arma::uword count = 0;
        arma::uvec indx = seq_int(window);
        // indx = indx-1;
        arma::uword number_of_blocks = floor(len/window);
        for ( arma::uword j = 0; j < number_of_blocks; ++j){
            RMS(i) = RMS(i) + arma::accu(arma::pow(poly_residuals(X.rows(indx), 
                                         order),2));
            count = count + 1;
            indx = indx + window;
        }
        
        RMS(i) = sqrt(RMS(i)/(count*window));

    }
    
    
    //take the logm of scales and RMS
    arma::vec log_scale(number_of_scales);
    arma::vec log_rms(number_of_scales);
    if (abs(scale_ratio-2)<.00001){
      log_scale = arma::log2(arma::conv_to<arma::vec>::from(scales));
      log_rms = arma::log2(RMS);
    }else{
      log_scale = log(arma::conv_to<arma::vec>::from(scales))/log(scale_ratio);
      log_rms = log(RMS)/log(scale_ratio);
    }

    //compute scaling coefficient
    arma::colvec alpha = lm_c(log_scale,log_rms);
    
    //create a list of output and return it
    if(verbose == 0){
        return List::create(Named("alpha") = alpha(1));
    }else{
        return List::create(Named("log_scales") = log_scale, Named("log_rms")=log_rms, Named("alpha") = alpha(1));
    }
}





//written by Aaron Likens (2022)
