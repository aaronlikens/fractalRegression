#' Multiscale Regression Plot
#' 
#' A plotting method for constructing scalewise regression plot
#' @param betas an object containing modeling results from multiscale regression
#' analysis. The object should be returned from the \code{mra} function of this
#' package.
#' @param order integer representing the detrending order used in the \code{mra} 
#' calculation. Default is 1.
#' @param ci a logical indicating whether confidence intervals should be 
#' computed using the \code{iaafft} function from this package. NOTE: with long 
#' time series (>> than N = 1,000), this can greatly reduce processing speed. 
#' Confidence intervals can be used for conventional significance testing of 
#' scale-wise correlation coefficients.
#' @param iterations integer that specifies the the number of surrogate time 
#' series to be generated for the purpose of confidence intervals. 
#' Default = 19. Larger number of surrogates will slow computational speed but
#' produce better confidence interval estimates.
#' @param return.ci logical indicating whether the confidence intervals 
#' should be returned
#' @export
mra.plot = function(betas, order = 1, ci = FALSE, iterations = NULL,
                    return.ci = FALSE){
  op = par(no.readonly = TRUE)
  par(mar = c(5,5,.5,1),
      cex.lab = 1.5,
      cex.axis = 1.5)
  
  
  if (ci){
    if (is.null(iterations)) {
      # for 95% confidence limit
      iterations = 19
    }
    
    # compute surrogate time series for x and y
    x.surr = iaafft(x, iterations)
    y.surr = iaafft(y, iterations)
    ci.betas = matrix(NA, nrow = iterations, ncol = length(betas$betas))
    cis = matrix(NA, nrow = 2, ncol = length(betas$scales))
    
    # perform dcca on each surrogagte
    for (i in 1:iterations){
      temp.x = as.vector(x.surr[,i])
      temp.y = as.vector(y.surr[,i])
      ci.betas[i,] = mra(x.surr[,i], y.surr[,i], order = order, 
                         scales = betas$scales)$betas
    }
    
    # compute 95% confidence intervals
    for (i in 1:length(betas$scales)){
      
      cis[,i] = quantile(ci.betas[,i], c(0.025, 0.9725))
    }
    ymin = min(min(ci.betas), min(betas$betas))
    ymax = max(max(ci.betas), max(betas$betas))
    plot(betas$scales, betas$betas, pch = 16, type = 'b', xlab = 's',
         ylab = expression(beta(s)), ylim = c(ymin, ymax))
    lines(betas$scales, cis[1,], col = 'red')
    lines(betas$scales, cis[2,], col = 'red')
  }
  par(op)
  
  if (return.ci){
    return(cis)
  }
  
}