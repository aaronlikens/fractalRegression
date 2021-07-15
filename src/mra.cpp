// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// these functions perform the DFA based regression introduced in:
// Kristoufek, L. (2015). 

// functions for doing the detrending covariance. Returns the sum of mean squared 
// DFA residuals and the covariance by scale
arma::vec detrend_cov(arma::vec x, arma::vec y, int m);
double detrendVar(arma::vec xr, int m);

//function for cumulative sums
arma::vec cumsum(arma::vec x);

//' Multiscale Regression Anlaysis (MRA)
//'
//' Fast function for computing multiscale regression analysis (MRA) on long time series. Combining DFA with ordinary least square regression, MRA
//' is a form of fractal regression that can be used to estimate asymmetric and multiscale regression coefficients between two variables. 
//'
//' @param x A real valued vector (i.e., time series data) to be analyzed. A key
//' difference between DCCA and MRA is that MRA produces asymmetric estiamtes. 
//' That is, x is assumed to be an independent variable and y is assumed to be 
//' a dependent variable. MRA should be used when one of the time series in 
//' question is usefully cast as the independent variable. That is, x is assumed
//' to effect change in y. If no such causal relationship is anticipated, use
//' DCCA instead.
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
//' @details Mathematical treatment of the MRA algorithm and its performance can be found in Kristoufek (2015) and Likens et al. (2019).
//'
//' Note that under conditions with linear and quadratic trends, Likens et al. (2019) found that there was a systematic positive bias in the \eqn{\beta} estimates for larger scales.
//' Using a polynomial detrending order of 2 or greater was shown to attenuate this bias. 
//'
//' The object returned from the mra() function will include the following:
//' \itemize{ 
//'  \item `scales` indicates the values of the scales used for estimates
//'  \item `betas` are the scale specific \eqn{\beta} estimates of the influence of x on y
//'  \item `r2` is the scale specific r-squared value of the model fit (i.e., variance in y accounted for by x at that scale)
//'  \item `t_observed` is the estimated t-statistic for a given \eqn{\beta} at a given scale. 
//' }
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
//' \dontrun{
//' # Here is a simple example for running MRA using a white noise and pink noise time series.
//' # For more detailed examples, see the vignette. 
//' 
//' noise <- rnorm(5000)
//' 
//' pink.noise <- fgn_sim(n = 5000, H = 0.9)
//'
//' scales <- ifultools::logScale(scale.min = 10, scale.max = 1250, scale.ratio = 1.1)
//' 
//' mra.out <- mra(x = noise, y = pink.noise, order = 1, scales = scales)
//' }
//'
//'
//[[Rcpp::export]]
List mra(arma::vec x, arma::vec y, int order, arma::ivec scales){
    int N = x.n_elem;
    double len = x.n_elem;
    int numberOfScales = scales.n_elem;
    arma::vec resid(numberOfScales);
    
    // create the profiles
    arma::vec X = cumsum(x-mean(x));
    arma::vec Y = cumsum(y-mean(y));
    
    // allocate vectors to fill with various scale by scale quantities
    // i.e., equations 4 - 8 in Kritoufek (2015).
     
    // scale by scale detrended variance in x
    arma::vec f2x(numberOfScales);

    // scale by scale detrended variance in y
    arma::vec f2y(numberOfScales);

    // scale by scale covariance
    arma::vec f2xy(numberOfScales);

    // error term scale by scale fluctuation
    arma::vec f2u(numberOfScales);

    //scale by scale error term, uhat_t(s)
    arma::vec ut(X.n_elem);

    //scale by scale R^2 values
    arma::vec r2(numberOfScales);

    //scale by scale regression coefficients
    arma::vec betas(numberOfScales);

    // standard error of scale by regression coefficient
    arma::vec se_betas(numberOfScales);

    // t-observed for the scale by scale regression coefficients
    arma::vec t_observed(numberOfScales);

    
    // Main loop for detrending and calculating the fluctuation functions
    for ( int i = 0; i < numberOfScales; i++){
        //choose window size
        int window = scales[i];
        //re-/initialize variables to hold input for a given scale
        int count = 0;
        arma::ivec indx = seq_len(window);
        indx = indx-1;
        arma::vec varCov(2);
        // std::fill(varCov.begin(),varCov.end(),0);
        int numberOfBlocks = floor(len/window);
        
        for ( int j = 0; j < numberOfBlocks; ++j){
            varCov = detrend_cov(X.subvec(indx(0), size(indx)), 
                                 Y.subvec(indx(0), size(indx)),
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
        betas[i] = f2xy[i]/f2x[i]; // scale-wise regression coefficient
        //betasYX = f2xy[i]/f2y[i];
        ut = Y - X*betas[i] - mean(Y-X*betas[i]);
        // //compute R^2 within each window and get an average
        indx = seq_len(window);
        indx = indx-1;
        //double beta = betas[i];
        for (int j = 0; j < numberOfBlocks; ++j){
            
            f2u[i] = f2u[i] + detrendVar(ut.subvec(indx(0), size(indx)), order);
            count = count + 1;
            indx = indx + window;
        }
        f2u[i] = f2u[i]/(N-window);
        se_betas[i] = sqrt(f2u[i]/(f2x[i]*(window-1)));
        t_observed[i] = betas[i]/se_betas[i];
        r2[i] = 1 - f2u[i]/f2y[i]; // scale wise 
        
                
    }
    
    return List::create(Named("scales") = scales, Named("betas") = betas,
                        Named("r2") = r2, Named("t_observed") = t_observed);
}


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

// function for doing the detrending. Returns the sum of squared residuals.
double detrendVar(arma::vec x, int m){
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
    
    // fit regression
    coefx = solve(t,x);
    
    // find residuals
    arma::colvec residx = x-t*coefx;
    
    //square and sum the residuals

    double f2x = 0;
    f2x = arma::accu(pow(residx,2))/(x.n_elem-1);
    
    return f2x;
}






//written by Aaron Likens (2019)
