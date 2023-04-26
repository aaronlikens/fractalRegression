// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// function for lagging variables
arma::mat LagN(arma::vec x, int lags);

// functions for detrending and computing scale-wise covariance
arma::vec PolyFit(arma::vec x, int order);
arma::mat Detrend_Var(arma::mat x, int order);

//' Multiscale Lagged Regression Anlaysis
//' Fast function for computing MLRA on long time series
//' @param x is a real valued vector of time series data
//' @param y is a real valued vector of time series data
//' @param order is an integer indicating the polynomial order used for 
//' detrending the local windows
//' @param scales integer vector of scales over which to compute correlation. 
//' Performance is best when scales are evenly spaced in log units. Choosing a
//' logarithm base between 1 and 2 may also improve performance of regression.
//' @param lags integer indicating the maximum number of lags to include
//' @param direction string indicating a positive ('p') or negative ('n') lag
//' @import Rcpp
//' @useDynLib fractalRegression
//' @export
//' 
//' @return The object returned from the dlcca() function is a list containing rho coefficients for each lag at each of the scales. 
//'
//[[Rcpp::export]]
arma::mat dlcca(arma::vec x, arma::vec y, int order, IntegerVector scales,
          int lags, char direction){
    
    // should be mat(N-l-1, l + 1), where l = number of lags
    arma::mat lx(x.n_elem - lags - 1, lags + 1, arma::fill::zeros);
    
    // should be mat(N-l-1, l + 1), where l = number of lags
    arma::mat ly(x.n_elem - lags - 1, lags + 1, arma::fill::zeros); 
    
    // should be mat(N-l-1, l + 2), where l = number of lags
    arma::mat lyx(x.n_elem - lags - 1, lags + 2, arma::fill::zeros); 
    
    unsigned int n_scales = scales.size();
    NumericVector resid(n_scales);
    arma::field<arma::mat> rmse(n_scales, 1);
    arma::mat betas(lags + 1, n_scales);

    arma::vec X = cumsum(x-mean(x));
    arma::vec Y = cumsum(y-mean(y));


    // negative
    if (direction == 'n'){
        lx = LagN(X, lags);
        ly = LagN(Y, lags);
        lyx.resize(ly.n_rows, ly.n_cols + 1);
        lyx.col(lyx.n_cols-1) = ly.col(0);
        for (unsigned int i = 0; i < lyx.n_cols-1; ++i){
            lyx.col(i) = lx.col(i);
        }
    // positive
    }else if (direction == 'p'){
        ly = LagN(Y, lags);
        lx = LagN(X, lags);
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
            
            varcov += Detrend_Var(lyx(arma::span(start,stop), arma::span::all), 
                                  order);
            start = start + scales[i];
            stop = stop + scales[i];
            count += 1;
            //Rcout << "Start = " << start <<". Stop = " << stop << "\n";
        }
        rmse(i,0) = varcov/(count*scales[i]);//check this for correctness later
        s = rmse(i,0);
        
        sxy = s(arma::span(0,s.n_rows-2), s.n_cols-1);
        sxx = s(arma::span(0,s.n_rows-2), arma::span(0,s.n_cols-2));
        // betas.col(i) = solve(sxx, sxy);
        // calculate scale-wise regression coefficients
        betas.col(i) = arma::inv(arma::diagmat(s))*s*arma::inv(arma::diagmat(s));
    }

    return betas;
}


arma::mat LagN(arma::vec x, int lags) {
    
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


arma::vec PolyFit(arma::vec x, int order){
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


arma::mat Detrend_Var(arma::mat X, int order){
    
    arma::mat d(X.n_rows, X.n_cols, arma::fill::zeros);
    
    for (unsigned int i = 0; i < X.n_cols; ++i){
        d.col(i) = PolyFit(X.col(i), order);
    }
    
    arma::mat varcov = cov(d);
    // arma::mat varcov = inv(d.t()*d);
    return varcov;
    
}
