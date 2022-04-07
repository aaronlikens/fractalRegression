#' Multifractional Brownian motion and multifractional Gaussian noise
#'
#' Simulate multifractional Brownian motion and multifractional Gaussian noise.
#'
#'
#'
#'
#' @param N The length of sample time series to simulate.
#' @param Ht The N by 1 vector of the time evolving H(t).
#'
#' @details This is an algorithm that simulates discrete time multifractional 
#' Brownian motion and multifractional Gaussian noise, which can useful for 
#' testing various functions within the `fractalRegression` package. H(t) 
#' should take on any values between 0 and 1. It is meant to capture time 
#' varying fractal properties. The example code given below shows a slow 
#' evolving Hurst exponent involving a sinusoidal change.  
#'
#' @return The object returned from the function includes:
#' \itemize{
#'  \item mBm: multifractional Brownian motion
#'  \item mGn: multifractional Gaussian noise
#' }
#'
#' @examples
#'
#' t <- 1:2048
#' Ht <- 0.5+0.5*(sin(0.0025*pi*t))
#' sim <- mBm_mGn(2048,Ht)
#'
#'
#'
#' @export


mBm_mGn <- function(N,Ht){


  numb1 <- 10
  numb2 <- 1000
  alpha <- 2
  N1 <- numb1*(numb2+N)

  mGnSum <- rep(0,numb1)
  mGnSumm <- rep(0,numb1*(numb2-1))
  mGn <- rep(0,N)

  R <- rnorm(N1)

  for (t in 1:N){
    for ( n in 1:numb1){
      mGnSum[n] <- (n^(Ht[t]-(1/alpha)))*R[1+(numb1*(numb2+t))-n]
    }
    mGnSum1 <- sum(mGnSum)
    for (nn in 1:(numb1*(numb2-1))){
      mGnSumm[nn] <- (((numb1+nn)^(Ht[t]-(1/alpha)))-nn^(Ht[t]-(1/alpha)))*R[1+(numb1*(numb2-1+t))-nn]
    }
    mGnSum2 <- sum(mGnSumm)
    mGn[t] <- ((numb1^(-Ht[t]))/gamma(Ht[t]-(1/alpha)+1))*(mGnSum1+mGnSum2)
  }
  mBm <- cumsum(mGn)
  out <- list(mBm,mGn)
  names(out) <- c("mBm", "mGn")
  return(out)
}
