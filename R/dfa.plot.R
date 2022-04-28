#' Detrended Fluctuation Plot
#' 
#' Plot method for monofractal detrended fluctuation analysis
#' 
#' 
#' 
#' @param x is an object returned from the \code{dfa} function of this package.
#' Plot parameters are chosen automatically, 
#' @export
#' 
dfa.plot = function(x){
  op = par(no.readonly = TRUE)
  par(mar = c(5,4.2,.5,1),
      cex.lab = 1.5,
      cex.axis = 1.5)
  m = lm(x$log_rms ~x$log_scales)
  plot(x$log_scales, x$log_rms, pch = 16, xlab = 'log(s)', ylab = 'logF(s)')
  abline(m, lwd = 2, col = 'red')
  eq.text = paste('logF(s) = ', format(coef(m)[1],digits = 2) ,' + ', 
                  format(coef(m)[2], digits = 2),
                  'log(s)' ,sep = '')
  r2 = format(summary(m)$r.square,4, digits = 3)
  r2.text = bquote(R^2*' = '~.(r2))
  text(min(x$log_scale)+1,max(x$log_rms)-.2, labels = eq.text)
  text(min(x$log_scale)+1, max(x$log_rms)-.4, labels = r2.text)
}