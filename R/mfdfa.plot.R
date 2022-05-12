#' Multifractal Spectrum Plot
#'
#'
#' Method for plotting various forms of the multifractal spectrum
#' 
#' 
#' @param mf an object containing elements related to the mutlifractal 
#' spectrum derived from Multifractal Detrended Fluctuation Analysis
#' @author Aaron D. Likens (2022)
#' @references Kantelhardt et al. (2002). Multifractal detrended fluctuation
#' analys of nonstationary time series. Physica A: Statistical Mechanics and 
#' its Applications, 87
#' @import colorRamps
#' @export
mfdfa.plot = function(mf){
  if (length(mf) != 7){
    cat('This does not appear to be a multifractal spectrum object.\n')
    return(NULL)
  }else{
    # gather current graphical device settings
    op = par(no.readonly = TRUE)
    
    # make some alterations to graphical device to adjust margins and font size
    par(mfrow = c(2,2),
        mar = c(5,4.2,.5,1),
        cex.lab = 1.5,
        cex.axis = 1.5)
    
    
    # do better color coding using a ramp function
    require(colorRamps)
    cols = rev(blue2red(length(mf$q)))
    
    # plot q-order fluctuation function
    # matplot(log10(mf$scales), log10(mf$fq),type = 'p',
    matplot(mf$log_scale, mf$log_fq,type = 'p',
            xlab = expression('logScale'), 
            ylab = expression('logF(q)'), 
            col = cols, tck = .03, pch = 16, cex = .75)
    
    # compute regressions and add regression lines to plot
    # regs = lm(log10(mf$fq) ~ log10(mf$scales))
    regs = lm(mf$log_fq ~ mf$log_scale)
    for (i in 1:length(mf$q)){
      abline(regs[[1]][,i], col = cols[i])
    }
    
    
    # plot the generalized Hurst exponents as a function of q
    plot(mf$q, mf$Hq, pch = 16,xlab = 'q', ylab = 'H(q)', col = cols,
         tck = .03, ylim = c(0.0, 1.2))
    
    # draw cross-hairs to to indicate the standard Hurst exponent derived from
    # monofractal detrended fluctuation anaysis
    abline(v = 2, col = 'red',lwd = 2)
    abline(h = mf$Hq[mf$q==2], col = 'red', lwd = 2)
    
    # plot the mass exponent as a function of q
    plot(mf$q, mf$tq, pch = 16, xlab = 'q', ylab = expression(tau(q)),
         col = cols, tck = .03)
    
    # plot the Legendre transformed multifractal spectrum
    plot(mf$h, mf$Dh, pch = 16, xlab = 'h', ylab = 'D(h)', 
         col = cols, tck = .03)
    
    # restore plot options to previous settings
    par(op)
  }
  
}