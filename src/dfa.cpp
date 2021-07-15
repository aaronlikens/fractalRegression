// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//detrending function that returns the sum of squared residuals
double polyDetrend(NumericVector yr, int m);

//simple linear regression
arma::colvec lmC(NumericVector xs, NumericVector ys);

//count the unique elements in a vector
int countUnique(NumericVector y);

//' Detrended Fluctuation Anlaysis
//' 
//' Fast function for computing detrended fluctuation analysis (DFA), a widely used method for estimating long-range temporal correlations in time series data. 
//' DFA is also a form of mono-fractal analysis that indicates the degree of self-similarity across temporal scales.
//' 
//' @param x A real valued vector (i.e., time series data) to be analyzed. 
//' @param order is an integer indicating the polynomial order used for 
//' detrending the local windows (e.g, 1 = linear, 2 = quadratic, etc.). There 
//' is not a pre-determined limit on the order of the polynomial order but the 
//' user should avoid using a large polynomial on small windows. This can result
//' in overfitting and non-meaningful estimates. 
//' @param verbose A logical that when = 1 indicates that the flucuation function including the log of all included scales as well as the log Rms should be 
//' returned as well as the \eqn{\alpha} or when = 0 only the estimated scaling exponent \eqn{\alpha} will be returned.
//' @param sc_min The minimum window size, specified in the number of data points (i.e., observations) to be included in the smallest window. 
//' @param sc_max The maximum window size, specified in the number of data points (i.e., observations) to be included in the largest window.
//' @param scale_ratio A scaling factor by which to create successive window sizes from `sc_min` to `sc_max`. 
//' This allows one to to maintain even spacing in logarithms while increasing
//' scale resolution.
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
//' The object returned from the function will include the following:
//' \itemize{ 
//'  \item If the value of verbose = 1, then a list object is returned that includes
//' the log of all included scales, the log root mean square error (RMS) per scale, and the overall \eqn{\alpha} estimate.
//'  \item If the value of verbose = 0, then only the estimated scaling exponent \eqn{\alpha} will be returned.
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
//' \dontrun{
//' 
//' noise <- rnorm(5000)
//' 
//' dfa.noise.out <- dfa(
//'     x = noise, 
//'     order = 1, 
//'     verbose = 1, 
//'     sc_min = 16, 
//'     sc_max = length(noise)/4, 
//'     scale_ratio = 2)
//' 
//' pink.noise <- fgn_sim(n = 5000, H = 0.9)
//' 
//' dfa.pink.out <- dfa(
//'     x = pink.noise, 
//'     order = 1, 
//'     verbose = 1, 
//'     sc_min = 16, 
//'     sc_max = length(pink.noise)/4, 
//'     scale_ratio = 2)
//' 
//' anticorr.noise <- fgn_sim(n = 5000, H = 0.25)
//' 
//' dfa.anticorr.out <- dfa(
//'     x = anticorr.noise, 
//'     order = 1, 
//'     verbose = 1, 
//'     sc_min = 16, 
//'     sc_max = length(anticorr.noise)/4, 
//'     scale_ratio = 2)
//' 
//' 
//' }
//' 
// [[Rcpp::export]]
List dfa(NumericVector x, int order, int verbose, double sc_min, 
         double sc_max, double scale_ratio ){
    
    double len = x.size();
    double scale_min = sc_min;
    double scale_max = sc_max;
    int numberOfScales = ceil(log(scale_max/scale_min)/log(scale_ratio));
    NumericVector resid(numberOfScales);
    NumericVector X = cumsum(x-mean(x));
    
    //create a logarithmiclly spaced vector of scales, contains duplicates
    NumericVector scales(numberOfScales);
    scales[0] = sc_min;
    for ( int i = 1; i < numberOfScales; ++i){
        scales[i] = scales[i-1]*scale_ratio;
    }
    //TODO (AARON): Find a way to sort the scales from lower to higher values
    scales = unique(floor(scales));
    
    // do the detrending and return the RMSE for each of the ith scales
    NumericVector RMS(numberOfScales);
    
    for ( int i = 0; i < numberOfScales; ++i){
        int window = scales[i];
        int count = 0;
        IntegerVector indx = seq_len(window);
        indx = indx-1;
        int numberOfBlocks = floor(len/window);
        for ( int j = 0; j < numberOfBlocks; ++j){
            RMS[i] = RMS[i] + polyDetrend(X[indx], order);
            count = count + 1;
            indx = indx + window;
        }
        
        RMS[i] = sqrt(RMS[i]/(count*window));
    }
    
    
    //take the logm of scales and RMS
    NumericVector logScale(numberOfScales);

    logScale = log(scales)/log(scale_ratio);
    NumericVector logRms(numberOfScales);
    logRms = log(RMS)/log(scale_ratio);

    //compute scaling coefficient
    arma::colvec alpha = lmC(logScale,logRms);
    
    //create a list of output and return it
    if(verbose == 0){
        return List::create(Named("alpha") = alpha(1));
    }else{
        return List::create(Named("logScales") = logScale, Named("logRms")=logRms, Named("alpha") = alpha(1));
    }
}

//detrending function that returns the sum of squared residuals
double polyDetrend(NumericVector yr, int m){
    int rows = yr.size();
    int cols = m + 1;
    arma::colvec coef(cols);
    
    arma::mat x(rows,cols);
    
    //convert data to an armadillo vector
    arma::colvec y(yr.begin(), yr.size(), false);
    
    //allocate memor for x and power of x vectors
    arma::colvec t1(rows);
    for ( int i = 0; i < rows; i++){
        t1(i) = i;
    }
    t1 = t1 - mean(t1);
    for ( int i = 0; i < cols; ++i){
        x.col(i) = arma::pow(t1,i);
    }
    
    coef = solve(x,y);
    
    NumericVector coefs(coef.begin(),coef.end());
    
    arma::colvec resid = y-x*coef;
    
    NumericVector Resid = NumericVector(resid.begin(),resid.end());
    
    //square and sum the residuals
    double ssResid = 0;
    for (int i = 0; i < rows; i ++){
        ssResid = ssResid + pow(Resid(i),2);
    }
    
    return ssResid;
}

//simple linear regression
arma::colvec lmC(NumericVector xs, NumericVector ys){
    NumericVector yr = ys;
    int n = yr.size();
    
    NumericVector Ones(n,1.0);
    NumericMatrix augmentedX(n,2);
    augmentedX(_,0) = Ones;
    augmentedX(_,1) = xs;
    
    arma::mat x(augmentedX.begin(), n, 2, false);
    arma::colvec y(yr.begin(), yr.size(), false);
    //arma::colvec w = 1/arma::linspace(1,n,n);
    //compute regression coefficients
    arma::colvec coef = solve(x,y);
    
    return coef;
}


// count the unique elements in a vector
int countUnique(NumericVector y) {
    NumericVector x = clone(y);
    std::sort(x.begin(), x.end());
    int count = 0;
    if (x[0] == x[1]){
        count += 1;
    }
    
    double diffx = 0.0;
    for ( int i = 1; i < x.size(); ++i){
        diffx = x[i] - x[i-1];
        if (diffx != 0){
            count +=1;
        }
    }
    
    return count;
}

//written by Aaron Likens (2019)
