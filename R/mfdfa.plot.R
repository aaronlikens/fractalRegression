#' Multifractal Spectrum Plot
#'
#'
#' Method for plotting various forms of the multifractal spectrum
#' 
#' 
#' @param mf an object containing elements related to the mutlifractal 
#' spectrum derived from Multifractal Detrended Fluctuation Analysis
#' @param do.surrogate logical indicating whether surrogation should be 
#' performed on the time series
#' @param nsurrogates integer indicating the number of surrogates to be 
#' constructed. Default is 19 for 95% confidence limit. Larger number of 
#' surrogates ore more precise but increase computational time.
#' @param return.ci logical indicating if confidence intervals derived from 
#' surrogate analysis should be returned. 
#' @author Aaron D. Likens (2022)
#' @references Kantelhardt et al. (2002). Multifractal detrended fluctuation
#' analys of nonstationary time series. Physica A: Statistical Mechanics and 
#' its Applications, 87
#' @importFrom colorRamps blue2red
#' @export
mfdfa.plot = function(mf, do.surrogate,  nsurrogates = 19, return.ci = FALSE){
  if (length(mf) != 11){
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
    # require(colorRamps)
    cols = rev(colorRamps::blue2red(length(mf$q)))
    
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
    
    
    # create surrogates and run multifractal analysis on each surrogate
    if (do.surrogate){
      x.surr = iaafft(mf$x, nsurrogates)
      x.surr.mfdfa = apply(x.surr, MARGIN = 2, FUN = function(x){
        out = mfdfa(x = x, q = mf$q, order = mf$order, scales = mf$scales,
                    scale_ratio = mf$scale_ratio)
      })
      
      # compute CI for each value of H(q)
      Hq.surr = lapply(x.surr.mfdfa, FUN =function(x){
        out = x$Hq
      })
      Hq.surr = Reduce(cbind, Hq.surr)
      
      Hq.ci = matrix(NA, nrow(Hq.surr), 2)
      for (i in 1:nrow(Hq.surr)){
        Hq.ci[i,] = as.vector(quantile(Hq.surr[i, ], c(0.025, 0.975)))
      }
      
      # compute CI for each value of tau(q)
      tau.surr = lapply(x.surr.mfdfa, FUN =function(x){
        out = x$Tau
      })
      
      tau.surr = Reduce(cbind, tau.surr)
      tau.ci = matrix(NA, nrow(tau.surr),ncol = 2)
      for (i in 1:nrow(tau.surr)){
        tau.ci[i,] = as.vector(quantile(tau.surr[i,], c(0.025, 0.975)))
      }
      
      # Compute CI for Holder exponent 
      h.surr = lapply(x.surr.mfdfa, FUN = function(x){
        out = x$h
      })
      h.surr = Reduce(cbind, h.surr)
      h.ci = matrix(NA, nrow(h.surr), ncol = 2)
      for (i in 1:nrow(h.surr)){
        h.ci[i, ] = as.vector(quantile(h.surr[i,], c(0.025, 0.975)))
      }
      
      # also get point estimate for h
      h.surr.ci.mean = rowMeans(h.surr)
      
      # compute CI for Holder dimension
      Dh.surr = lapply(x.surr.mfdfa, FUN = function(x){
        out = x$Dh
      })
      Dh.surr = Reduce(cbind, Dh.surr)
      Dh.ci = matrix(NA, nrow(Dh.surr), ncol = 2)
      for (i in 1:nrow(Dh.surr)){
        Dh.ci[i, ] = as.vector(quantile(Dh.surr[i, ], c(0.025, 0.975)))
      }
      
      # also get point estimates for Dh
      Dh.surr.ci.mean = rowMeans(Dh.surr)
      
      # plot H(q) with confidence intervals as a line
      hq.ymax = max(max(Hq.ci[,2]), max(mf$Hq))
      hq.ymin = min(min(Hq.ci[,1]), min(mf$Hq))
      plot(mf$q, mf$Hq, pch = 16,xlab = 'q', ylab = 'H(q)', col = cols,
           tck = .03, ylim = c(hq.ymin, hq.ymax))
      
      lines(mf$q, Hq.ci[,1], col = 'black', lty = 2)
      lines(mf$q, Hq.ci[,2], col = 'black', lty = 2)
      abline(v = 2, col = 'red',lwd = 2)
      abline(h = mf$Hq[mf$q==2], col = 'red', lwd = 2)
      
      
      
      # plot the mass exponent as a function of q
      tau.ymax = max(max(tau.ci[,2]), max(mf$Tau))
      tau.ymin = min(min(tau.ci[,1]), min(mf$Tau))
      plot(mf$q, mf$Tau, pch = 16, xlab = 'q', ylab = expression(tau(q)),
           col = cols, tck = .03, ylim = c(tau.ymin, tau.ymax))
      lines(mf$q, tau.ci[,1], col = 'black', lty = 2)
      lines(mf$q, tau.ci[,2], col = 'black', lty = 2)
      
      # plot the Legendre transformed multifractal spectrum
      hmin = min(min(h.ci[,1]), min(mf$h))
      hmax = max(max(h.ci[,2]), max(mf$h))
      Dhmin = min(min(Dh.ci[,1]), min(mf$Dh))
      Dhmax = max(max(Dh.ci[,2]), max(mf$Dh))
      plot(mf$h, mf$Dh, pch = 16, xlab = 'h', ylab = 'D(h)', 
           col = cols, tck = .03, xlim = c(hmin, hmax), ylim = c(Dhmin, Dhmax))
      
      # add confidence intervals for dh/Dh pairs
      points(h.surr.ci.mean, Dh.surr.ci.mean, pch = 16, cex = 0.5, 
             col = cols)
      suppressWarnings(
      arrows(x0 = h.ci[,1], y0 = Dh.surr.ci.mean,
             x1 = h.ci[,2], y1 = Dh.surr.ci.mean,
             code = 3, angle = 90, length = 0.05,
             col = cols)
      )
      suppressWarnings(
      arrows(x0 = h.surr.ci.mean, y0 = Dh.ci[,1],
             x1 = h.surr.ci.mean, y1 = Dh.ci[,2],
             code = 3, angle = 90, length = 0.05, 
             col = cols)
      )
    }else{
      plot(mf$q, mf$Hq, pch = 16,xlab = 'q', ylab = 'H(q)', col = cols,
           tck = .03)
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
    }
    
    

    
    # restore plot options to previous settings
    par(op)
  }
  
}