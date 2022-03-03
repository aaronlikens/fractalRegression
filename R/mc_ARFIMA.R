#' Mixed-correlated ARFIMA processes
#'
#' Simulate various types of correlated noise processes.
#'
#' @param process specifies the type of correlated noise process to simulate
#' and includes 'Noise_rho', 'ARFIMA_rho','ARFIMA_AR','AR_rho',
#' 'Mixed_ARFIMA_ARFIMA','Mixed_ARFIMA_AR',and 'Mixed_ARFIMA_noise'.
#' @param n is a numeric value specifying the length of the time-series.
#' @param rho specifies the strength of the correlation with values -1 - 1.
#' @param d1 is a numeric fractional difference parameter for x specifying long term 
#' memory.
#' @param d2 is a numeric fractional difference parameter for x specifying long term 
#' memory.
#' @param d3 is a numeric fractional difference parameter for y specifying long term 
#' memory.
#' @param d4 is a numeric fractional difference parameter for y specifying long term 
#' memory.
#' @param alpha
#' @param beta
#' @param delta
#' @param gamma
#' @param theta
#' @param theta1
#' @param theta2
#'
#'
#'
#' @details This function includes multiple options simulating various
#' types of correlated noise processes including mixed-correlated ARFIMA
#' processes with power-law cross-correlations, These functions were originally
#' written by Ladislav Kristoufek and posted on his website. They go with
#' the paper presented in Kristoufek (2013). The 'process' argument specifies
#' the type of noise to be generated.
#'  \itemize{
#'   \item 'Noise_rho' - Generates two correlated noise series and requires
#'   arguments: n, rho.
#'   \item 'ARFIMA_rho' - Generates two ARFIMA processes with correlated
#'   innovations and requires arguments: n, d1, d2, rho.
#'   \item 'ARFIMA_AR'  - Generates ARFIMA and AR(1) processes with correlated
#'   innovations and requires arguments: n, d1, theta, rho.
#'   \item 'AR_rho' - Generates two AR(1) processes with correlated innovations
#'   and requires arguments: n, theta1, theta2, rho.
#'   \item 'Mixed_ARFIMA_ARFIMA' - Generates MC-ARFIMA process with long-range
#'   correlation and long-range cross-correlation (Kristoufec, 2013 Model 1) and
#'    requires arguments: alpha, beta, gamma, delta, n, d1, d2, d3, d4, rho.
#'   \item 'Mixed_ARFIMA_AR' - Generates MC-ARFIMA process with long-range
#'   correlation and short-range cross-correlation (Kristoufec, 2013 Model 2)
#'   and requires arguments: alpha, beta, gamma, delta, n, d1, d2, theta, rho.
#'   \item 'Mixed_ARFIMA_noise' - Generates MC-ARFIMA process with long-range
#'   correlation and simple correlation (Kristoufec, 2013 Model 3) and requires
#'   arguments: alpha, beta, gamma, delta, n, d1, d2, rho.
#'  }
#'
#' @return
#' The object returned is a matrix of length n with a time series (x,y)
#' in column 1 and 2.
#'
#' @references
#' Kristoufek, L. (2013). Mixed-correlated ARFIMA processes for power-law
#' cross-correlations. Physica A: Statistical Mechanics and its Applications,
#' 392(24), 6484-6493.
#'
#'
#'
#'
#' @examples
#' set.seed(987345757)
#'
#' sim1 <- mc_ARFIMA(process='Mixed_ARFIMA_ARFIMA', alpha = 0.2,
#' beta = 1, gamma = 1, delta = 0.2, n = 10000, d1 = 0.4, d2 = 0.3,
#' d3 = 0.3, d4=0.4, rho=0.9)
#'
#' plot(sim1[,1],type='l', ylab= "Signal Amplitude", xlab='Time',
#' main = "MC-ARFIMA with LRC and LRCC")
#'
#' lines(sim1[,2], col='blue')
#'


mc_ARFIMA <- function(process,n, rho, d1=NULL,d2=NULL,d3=NULL,d4=NULL,alpha=NULL,beta=NULL,delta=NULL,gamma=NULL,theta=NULL,theta1=NULL,theta2=NULL){
  require(fracdiff)
  # Correlated Noise
  if(process=='Noise_rho'){
  e1<-rnorm(n)
  e<-rnorm(n)
  e2<-rho*e1+sqrt(1-rho^2)*e
  return(matrix(c(e1,e2),ncol=2))
  #Two ARFIMA processes with correlated innovations
  } else if(process=="ARFIMA_rho") {
    e1<-rnorm(n)
    e<-rnorm(n)
    e2<-rho*e1+sqrt(1-rho^2)*e
    x<-fracdiff.sim(n,d=d1,innov=e1)$series
    y<-fracdiff.sim(n,d=d2,innov=e2)$series
    return(matrix(c(x,y),ncol=2))
  #ARFIMA and AR(1) processes with correlated innovations
  } else if(process=="ARFIMA_AR") {
    e1<-rnorm(n)
    e<-rnorm(n)
    e2<-rho*e1+sqrt(1-rho^2)*e
    x<-fracdiff.sim(n,d=d1,innov=e1)$series
    y<-arima.sim(list(order=c(1,0,0),ar=theta),n,innov=e2)
    return(matrix(c(x,y),ncol=2))
#Two AR(1) processes with correlated innovations
  } else if(process=="AR_rho") {
    e1<-rnorm(n)
    e<-rnorm(n)
    e2<-rho*e1+sqrt(1-rho^2)*e
    x<-arima.sim(list(order=c(1,0,0),ar=theta1),n,innov=e1)
    y<-arima.sim(list(order=c(1,0,0),ar=theta2),n,innov=e2)
    return(matrix(c(x,y),ncol=2))
#MC-ARFIMA process with LRC and LRCC (Model 1)
  } else if(process=="Mixed_ARFIMA_ARFIMA") {
    er1<-rnorm(n)
    er<-rnorm(n)
    er2<-rho*er1+sqrt(1-rho^2)*er
    x1<-fracdiff.sim(n,d=d1,innov=er1)$series
    y1<-fracdiff.sim(n,d=d2,innov=er2)$series
    e<- matrix(c(x1,y1),ncol=2)
    e2<-e[,1]
    e3<-e[,2]
    x<-alpha*fracdiff.sim(n,d=d1)$series+beta*e2
    y<-gamma*e3+delta*fracdiff.sim(n,d=d4)$series
    return(matrix(c(x,y),ncol=2))
    #MC-ARFIMA process with LRC and SRCC (Model 2)
  } else if(process=="Mixed_ARFIMA_AR") {
    er1<-rnorm(n)
    er<-rnorm(n)
    er2<-rho*er1+sqrt(1-rho^2)*er
    x1<-arima.sim(list(order=c(1,0,0),ar=theta1),n,innov=er1)
    y1<-arima.sim(list(order=c(1,0,0),ar=theta2),n,innov=er2)
    e<-matrix(c(x1,y1),ncol=2)
    e2<-e[,1]
    e3<-e[,2]
    x<-alpha*fracdiff.sim(n,d=d1)$series+beta*e2
    y<-gamma*e3+delta*fracdiff.sim(n,d=d2)$series
    return(matrix(c(x,y),ncol=2))
  #MC-ARFIMA process with LRC and simple correlation (Model 3)
    } else {
      er1<-rnorm(n)
      er<-rnorm(n)
      er2<-rho*er1+sqrt(1-rho^2)*er
      e<-matrix(c(er1,er2),ncol=2)
      e2<-e[,1]
      e3<-e[,2]
      x<-alpha*fracdiff.sim(n,d=d1)$series+beta*e2
      y<-gamma*e3+delta*fracdiff.sim(n,d=d2)$series
      return(matrix(c(x,y),ncol=2))
    }

  }

