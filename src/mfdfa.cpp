// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//helper function prototypes
NumericVector poly_residuals(NumericVector yr, int m);  //use polynomial detrending and return residuals
NumericVector logscale(double smin, double smax, double sratio); //create logarithmically spaced series
NumericVector lm_c(NumericVector xs, NumericVector yr); // simple first order regression
NumericVector linspace(double from, double to, int length_out); //create linearly increasing series
NumericVector colmeans(NumericMatrix X); //find the mean for each column
int which_value(NumericVector x, double value); //find the first index where a value == to a specified number

//' Multifractal Detrended Fluctuation Analysis
//'
//' Fast function for computing multifractal detrended fluctuation analysis (MF-DFA), a widely used method for estimating the family of long-range temporal correlations or scaling exponents in time series data. 
//' MF-DFA is also a form of multifractal analysis that indicates the degree of interaction across temporal scales.
//' 
//' @param x A real valued vector (i.e., time series data) to be analyzed. 
//' @param q A real valued vector indicating the statistical moments (q) to use 
//' in the analysis.
//' @param order is an integer indicating the polynomial order used for 
//' detrending the local windows (e.g, 1 = linear, 2 = quadratic, etc.). There 
//' is not pre-determined limit on the order of the polynomial order but the 
//' user should avoid using a large polynomial on small windows. This can result
//' in overfitting and non-meaningful estimates. 
//' @param scale_min An integer indicating the minimum window size, specified in the number of data points (i.e., observations) to be included in the smallest window. 
//' @param scale_max An integer indicating the maximum window size, specified in the number of data points (i.e., observations) to be included in the largest window. indicating largest scale to resolve
//' @param scale_ratio A scaling factor by which to create successive window sizes from `scale_min` to `scale_max`. 
//' This allows one to to maintain even spacing in logarithms while increasing
//' scale resolution.
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
//' The output of the algorithm is a list that includes:
//' \itemize{ 
//'  \item log2scale: The log2 scales used for the analysis
//'  \item log2Fq: The log2 of the fluctuation functions for each scale and q 
//'  \item Hq: The q-order Hurst exponent (generalized Hurst exponent)
//'  \item Tau: The q-order mass exponent
//'  \item q: The q-order statistical moments
//'  \item h: The q-order singularity exponent
//'  \item Dh: The dimension of the q-order singularity exponent
//'}
//' While it is common to use only linear detrending with DFA and MF-DFA, it is important to inspect the trends in the data to determine
//' if it would be more appropriate to use a higher order polynomial for detrending, and/or compare the DFA and MF-DFA output for different polynomial orders (see Ihlen, 2012; Kantelhardt et al., 2001).
//' 
//' General recommendations for choosing the min and max scale are a scale_min = 10 and scale_max = (N/4), where N is the number of observations.
//' See Eke et al. (2002), Gulich and Zunino (2014), Ihlen (2012), and  for additional considerations and information on choosing the correct parameters. 
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
//' \dontrun{
//' 
//' noise <- rnorm(5000)
//' 
//' mf.dfa.white.out <- mfdfa(
//'     x = noise, q = c(-5:5), 
//'     order = 1, 
//'     scale_min = 16, 
//'     scale_max = length(noise)/4, 
//'     scale_ratio = 2) 
//'  
//' pink.noise <- fgn_sim(n = 5000, H = 0.9)
//' 
//' mf.dfa.pink.out <- mfdfa(
//'     x = pink.noise, 
//'     q = c(-5:5), 
//'     order = 1, 
//'     scale_min = 16, 
//'     scale_max = length(pink.noise)/4, 
//'     scale_ratio = 2)
//'
//' }
//' 
// [[Rcpp::export]]
List mfdfa(NumericVector x, NumericVector q, int order, double scale_min, 
           double scale_max, double scale_ratio){
    try{
        //double scale_ratio = 1.1;
        double len = x.size();  // get size of time series
        //double scale_min = 256;  //minimum scale to resolve
        int numberOfScales = ceil(log(scale_max/scale_min)/log(scale_ratio));//determine how many scales to use
        NumericVector X = cumsum(x-mean(x));  //take the cumulative sum of the data
        
        //create vectors of scales and q values
        NumericVector scales = logscale(scale_min, scale_max, scale_ratio);
        int qlength = q.size();
        int q0indx = which_value(q, 1e-13);
        q[q0indx] = 0.0;
        
        //do the detrending and return the RMSE for each of the ith scales
        NumericMatrix fq(numberOfScales,qlength);
        NumericMatrix log2Fq(numberOfScales,qlength);
        NumericMatrix qRMS(numberOfScales, qlength);
        
        for ( int ns = 0; ns < numberOfScales; ++ns ){
            int window = scales[ns];
            IntegerVector indx(window);
            for (int j = 0; j < window; ++j){
                indx[j] = j;
            }
            int numberOfBlocks = floor(len/window);
            NumericVector rms(numberOfBlocks);
            NumericMatrix qrms(numberOfBlocks,qlength);
            
            for ( int v = 0; v < numberOfBlocks; ++v ){
                rms[v] = mean(pow(poly_residuals(X[indx],order),2));
                rms[v] = sqrt(rms[v]);
                indx = indx + window;
                for ( int nq = 0; nq < qlength; ++nq ){
                    qrms(v,nq) = pow(rms[v],q[nq]);
                }
            }
            
            qRMS(ns, _) = colmeans(qrms); //compute the q-order statistics
            
            for ( int nq = 0; nq < qlength; ++nq){
                fq(ns,_) = pow(qRMS(ns,_),(1/q[nq]));
                log2Fq(ns,nq) = log2(fq(ns,nq));
            }
            
            log2Fq(ns,q0indx) = (log2Fq(ns,q0indx-1) + log2Fq(ns,q0indx+1))/2;
        }
        
        //take the log2 of scales
        NumericVector log2Scale(numberOfScales);
        for ( int i = 0; i < numberOfScales; ++i ){
            log2Scale[i] = log2(scales[i]);
        }
        
        //compute various fractal scaling exponents (tau, h, Dh)
        NumericVector Hq(q.size());
        
        for ( int nq = 0; nq < q.size(); ++nq ){
            NumericVector temp = log2Fq(_,nq);
            NumericVector p = lm_c(log2Scale,temp);
            Hq[nq]=p[1];
        }
        
        NumericVector tau(Hq.size());
        tau = Hq*q-1;
        NumericVector hh = diff(tau)/(q[1]-q[0]); 
        NumericVector Dh(q.size()-1);
        
        for ( int i = 0; i < q.size(); ++i ){
            Dh[i] = q[i]*hh[i]-tau[i];
        }
        
        NumericVector h = hh-1;
        
        return List::create(Named("log2Scale") = log2Scale, Named("log2Fq")=log2Fq, Named("Hq") = Hq,
        Named("Tau") = tau, Named("q") = q, Named("h")= h,
        Named("Dh")=Dh);
        
    } catch( std::exception &ex ) {    	// or use END_RCPP macro
    forward_exception_to_r( ex );
    } catch(...) { 
        ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue; // -Wall
}

//function to compute the residuals
NumericVector poly_residuals(NumericVector yr, int m){
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
    
    for ( int i = 0; i < cols; ++i){
        x.col(i) = arma::pow(t1,i);
    }
    
    coef = solve(x,y);
    
    NumericVector coefs(coef.begin(),coef.end());
    
    arma::colvec resid = y-x*coef;
    
    NumericVector Resid = NumericVector(resid.begin(),resid.end());
    
    return Resid;
}

//function for creating an exponentially increasing vector
NumericVector logscale(double smin, double smax, double sratio){
    double NumberOfScales = ceil(log(smax/smin)/log(sratio));
    NumericVector Scales(NumberOfScales);
    Scales[0] = smin;
    for ( int i = 1; i< NumberOfScales; i++){
        Scales[i] = Scales[i-1]*sratio;
    }
    return Scales;
}

//compute simple linear slope and intercept
NumericVector lm_c(NumericVector xs, NumericVector yr){
    NumericVector coefs;
    int n = yr.size();
    
    NumericVector Ones(n,1.0);
    NumericMatrix augmentedX(n,2);
    augmentedX(_,0) = Ones;
    augmentedX(_,1) = xs;
    
    arma::mat x(augmentedX.begin(), n, 2, false);
    arma::colvec y(yr.begin(), yr.size(), false);
    
    //compute regression coefficients
    arma::colvec coef = solve(x,y);
    
    coefs = NumericVector(coef.begin(),coef.end());
    return coefs;
}

//create a linear sequence
NumericVector linspace(double from, double to, int length_out){
    double len = length_out;
    double increment = (to-from)/(len-1);
    NumericVector vec(length_out);
    vec[0]=from;
    for ( int i = 0; i < length_out; i++){
        vec[i+1] = vec[i] + increment;
    }
    return vec;
}

//find the column means of a matrix
// [[Rcpp::export]]
NumericVector colmeans(NumericMatrix X){
    int cols = X.ncol();
    NumericVector means(cols);
    for ( int i = 0; i < cols; i++){
        means[i] = mean(X(_,i));
    }
    return means;
}

//find a specific value. If exported and used in R, remember to add one to the result to
//account for the indexing differences in R vs. c++
int which_value(NumericVector x, double value){
    LogicalVector res = x <= value;
    int indx = 0;
    for ( int i = 0; i < x.size(); i++){
        if(res[i] == TRUE){
            indx = i;
        }
    }
    return indx;
}

//written by Aaron Likens (2019)
