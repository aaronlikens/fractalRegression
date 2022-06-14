// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "lmc.h"
using namespace Rcpp;
    

//' Detrended Cross-Correlation Analysis
//'
//' Fast function for computing detrended cross-correlation analysis (DCCA) on long time series, which is a 
//' bivariate extension of detrended fluctuation analysis (DFA).
//'
//' @param x A real valued vector (i.e., time series data) to be analyzed.
//' @param y A real valued vector (i.e., time series data) to be analyzed.
//' @param order is an integer indicating the polynomial order used for 
//' detrending the local windows (e.g, 1 = linear, 2 = quadratic, etc.). There 
//' is not a pre-determined limit on the order of the polynomial order but the 
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
//' @import Rcpp
//' @useDynLib fractalRegression
//' @export
//'
//' @details Details of the algorithm are specified in Podobnik and Stanley (2008) and in Zebende (2011).
//' In general, the output of the algorithm are estimates of \eqn{\rho}DCCA, which range from -1 to 1 and can generally be interpreted as follows: 
//' \itemize{
//'  \item \eqn{\rho DCCA = -1.0 ->} perfect anti-cross-correlation
//'  \item \eqn{\rho DCCA =  0.0 ->} no cross-correlation
//'  \item \eqn{\rho DCCA =  1.0 ->} perfect cross-correlation
//' } 
//'
//' @return The object returned from the function is a list including the following:
//' \itemize{ 
//'  \item \code{scales} indicates the values of the scales used for estimates \eqn{\rho}DCCA
//'  \item \code{rho} includes the scale-wise estimates of \eqn{\rho}DCCA
//' }
//'
//' @references
//' 
//' Podobnik, B., & Stanley, H. E. (2008). Detrended cross-correlation analysis: a new method for analyzing two nonstationary time series. Physical review letters, 100(8), 084102.
//'
//' Zebende, G. F. (2011). DCCA cross-correlation coefficient: Quantifying level of cross-correlation. Physica A: Statistical Mechanics and its Applications, 390(4), 614-618.
//'
//'//'
//' @examples
//' 
//'
//' 
//' # Here is a simple example for running DCCA using a white noise and pink noise time series.
//' # For more detailed examples, see the vignette. 
//' 
//' noise <- rnorm(5000)
//' 
//' pink.noise <- fgn_sim(n = 5000, H = 0.9)
//'
//' scales <- ifultools::logScale(scale.min = 10, scale.max = 1250, scale.ratio = 1.1)
//' 
//' dcca.out <- dcca(noise, pink.noise, order = 1, scales = scales)
//' 
//' 
//'
// [[Rcpp::export]]
List dcca(arma::vec x, arma::vec y, int order, arma::ivec scales){
    int N = x.n_elem;
    double len = x.n_elem;
    int numberOfScales = scales.n_elem;
    arma::vec resid(numberOfScales);
    arma::uword start = 0;
    arma::uword stop = 0;
    
    // create the profiles
    arma::vec X = cumsum(x-mean(x));
    arma::vec Y = cumsum(y-mean(y));
    
    // allocate vectors to fill with various scale by scale quantities
    // i.e., equations 4 - 8 in Kritoufek (2015).
     
    // scale by scale detrended variance in x
    arma::vec f2x(numberOfScales);
    //std::fill(f2x.begin(),f2x.end(),0);
    
    // scale by scale detrended variance in y
    arma::vec f2y(numberOfScales);
    //std::fill(f2y.begin(),f2y.end(),0);
    
    // scale by scale covariance
    arma::vec f2xy(numberOfScales);
    //std::fill(f2xy.begin(),f2xy.end(),0);
    
    //scale by scale regression coefficients
    arma::vec rho(numberOfScales);
    //std::fill(rho.begin(),rho.end(),0.0);
    
    // Main loop for detrending and calculating the fluctuation functions
    for ( int i = 0; i < numberOfScales; i++){
        //choose window size
        int window = scales[i];
        //re-/initialize variables to hold input for a given scale
        int count = 0;
        // IntegerVector indx = seq_len(window);
        arma::ivec indx = seq_len(window);
        indx = indx-1;
        arma::vec varCov(2);
        //std::fill(varCov.begin(),varCov.end(),0);
        int numberOfBlocks = floor(len/window);
        
        for ( int j = 0; j < numberOfBlocks; ++j){
            // varCov = detrend_cov(X[indx],Y[indx],order);
            // varCov = detrend_cov(X.subvec(indx(0), size(indx)), Y.subvec(indx(0),
            //                        size(indx)), order);
            start = indx(0);
            stop = indx(indx.n_elem-1);
            varCov = detrend_cov(X.subvec(start, stop),
                               Y.subvec(start, stop),
                               order);
            f2x[i] = f2x[i] + varCov[0];
            f2y[i] = f2y[i] + varCov[1];
            f2xy[i] = f2xy[i] + varCov[2];
            count = count + 1;
            indx = indx + window;
        }

        f2x[i] = f2x[i]/(N-window);
        f2y[i] = f2y[i]/(N-window);
        f2xy[i] = f2xy[i]/(N-window);
        
        // scale-wise correlation coefficient
        rho[i] = f2xy[i]/(sqrt(f2x[i])*sqrt(f2y[i])); 

        
                
    }
    
    return List::create(Named("x") = x,
                        Named("y") = y,
                        Named("scales") = scales,
                        Named("rho") = rho);
}





//written by Aaron Likens (2019)
