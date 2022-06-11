// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "lmc.h"
using namespace Rcpp;


//' Multifractal Detrended Fluctuation Analysis
//'
//' Fast function for computing multifractal detrended fluctuation analysis 
//' (MF-DFA), a widely used method for estimating the family of long-range 
//' temporal correlations or scaling exponents in time series data. 
//' MF-DFA is also a form of multifractal analysis that indicates the degree 
//' of interaction across temporal scales.
//' 
//' @param x A real valued vector (i.e., time series data) to be analyzed. 
//' @param q A real valued vector indicating the statistical moments (q) to use 
//' in the analysis. q must span negative and positive values e.g., -3:3, 
//' otherwise and error may be produced. 
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
//'  \item \code{log_scale} The log2 scales used for the analysis
//'  \item \code{log_fq} The log2 of the fluctuation functions for each scale and q 
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
List mfdfa(arma::vec x, arma::vec q, int order, arma::uvec scales){
    try{
        double len = x.size();  // get size of time series
        unsigned int numberOfScales = scales.n_elem;//determine how many scales to use
        arma::vec X = cumsum(x-mean(x));  //take the cumulative sum of the data
        
        //create vectors of scales and q values
        unsigned int qlength = q.n_elem;
        arma::uvec q0indx = arma::find(abs(q) <= 1e-8, 1); // may be incompatible
        bool q_contains_zero = arma::numel(q0indx) > 0;
        
        if (q_contains_zero) {
          q(q0indx(0)) = 0.0;  
        }
        
        
        //do the detrending and return the RMSE for each of the ith scales
        arma::mat fq(numberOfScales,qlength);
        arma::mat log_fq(numberOfScales,qlength);
        arma::mat qRMS(numberOfScales, qlength);
        
        for (unsigned int ns = 0; ns < numberOfScales; ++ns ){
            unsigned int window = scales(ns);
            arma::uvec indx(window);
            for (unsigned int j = 0; j < window; ++j){
                indx(j) = j;
            }
            unsigned int numberOfBlocks = floor(len/window);
            arma::vec rms(numberOfBlocks);
            arma::mat qrms(numberOfBlocks,qlength);
            arma::vec resid(indx.n_elem); 
            arma::vec resid_sq(indx.n_elem); 
           
            for (unsigned int v = 0; v < numberOfBlocks; ++v ){
              
                arma::vec temp = X.rows(indx); 
                resid = poly_residuals(temp,order);
                resid_sq = arma::pow(resid,2);
                rms[v] = arma::mean(resid_sq);
                rms[v] = sqrt(rms[v]);
                indx = indx + window;
                for (unsigned int nq = 0; nq < qlength; ++nq ){
                    qrms(v,nq) = pow(rms(v),q(nq));
                }
                
            }
            
            qRMS.row(ns) = mean(qrms,0); //compute the q-order statistics
            
            for (unsigned int nq = 0; nq < qlength; ++nq){
                fq.row(ns) = arma::pow(qRMS.row(ns),(1/q[nq]));
                log_fq(ns,nq) = log2(fq(ns,nq));

            }
            // TODO: address issue where q does not contain negative values
            if (q_contains_zero){
              log_fq(ns, q0indx(0)) = (log_fq(ns, q0indx(0)-1) + log_fq(ns, q0indx(0)+1))/2;  
            }
            
        }
        
        
        //take the log2 of scales
        arma::vec log_scale(numberOfScales);
        for ( unsigned int i = 0; i < numberOfScales; ++i ){
            log_scale(i) = log2(scales(i));
        }
        
        //compute various fractal scaling exponents (tau, h, Dh)
        arma::vec Hq(q.n_elem);
        for ( unsigned int nq = 0; nq < q.n_elem; ++nq ){
            arma::vec temp = log_fq.col(nq);
            arma::vec p = lm_c(log_scale,temp);
            Hq(nq)=p(1);
        }

        
        arma::vec tau(Hq.n_elem);
        tau = Hq%q-1;
        arma::vec tau_diff = arma::diff(tau);
        double q_increment = q(1) - q(0);
        arma::vec hh = tau_diff/q_increment; 
        arma::vec Dh(q.n_elem-1);
        
        for (unsigned int i = 0; i < q.n_elem-1; ++i ){
            Dh(i) = q(i)*hh(i)-tau(i);
        }
        
        arma::vec h = hh;
        
        return List::create(Named("log_scale") = log_scale, Named("log_fq")=log_fq, Named("Hq") = Hq,
        Named("Tau") = tau, Named("q") = q, Named("h")= h,
        Named("Dh")=Dh);
        
    } catch( std::exception &ex ) {    	// or use END_RCPP macro
    forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "c++ exception (unknown reason)" );
    }
    return R_NilValue; // -Wall
}

//written by Aaron Likens (2022)
