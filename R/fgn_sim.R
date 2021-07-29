# this software was taken from the now defunct R package
# fArma.  

#' Simulate fractional Gaussian Noise.
#' @param n integer indicating length of desired series
#' @param H Hurst exponent ranges between 0 and 1
#' @importFrom stats fft rnorm
#' @return A numeric vector of length n.
#' @export
fgn_sim <- 
    function(n = 1000, H = 0.7)
    {   
        # A function implemented by Diethelm Wuertz
        
        # Description:
        #   Creates a series of fractional Gaussian Noise
        
        # Arguments:
        #   n - length of desired series.
        #   H - the Hurst parameter.
        #   sigma - the standard deviation of the innovations used.
        
        # Details:
        #   FGN time series simulation. The FGN sequences use
        #   functions going back to Taqqu et al., Paxson and
        #   Beran (with Maechler's modifications from StatLib).
        
        # FUNCTION:
        
        # Settings:
        mean = 0
        std = 1
        ans = NA
        # Generate Sequence:
        z = rnorm(2*n)
        zr = z[c(1:n)]
        zi = z[c((n+1):(2*n))]
        zic = -zi
        zi[1] = 0
        zr[1] = zr[1]*sqrt(2)
        zi[n] = 0
        zr[n] = zr[n]*sqrt(2)
        zr = c(zr[c(1:n)], zr[c((n-1):2)])
        zi = c(zi[c(1:n)], zic[c((n-1):2)])
        z = complex(real = zr,imaginary = zi)
        
        # .gkFGN0:
        k = 0:(n-1)
        gammak = (abs(k-1)**(2*H)-2*abs(k)**(2*H)+abs(k+1)**(2*H))/2
        ind = c(0:(n - 2), (n - 1), (n - 2):1)
        .gkFGN0 = fft(c(gammak[ind+1]), inverse = TRUE)
        gksqrt = Re(.gkFGN0)
        if (all(gksqrt > 0)) {
            gksqrt = sqrt(gksqrt)
            z = z*gksqrt
            z = fft(z, inverse = TRUE)
            z = 0.5*(n-1)**(-0.5)*z
            z = Re(z[c(1:n)])
        } else {
            gksqrt = 0*gksqrt
            stop("Re(gk)-vector not positive")
        }
        
        # Standardize:
        # (z-mean(z))/sqrt(var(z))
        ans = std*drop(z) + mean
        return(ans)
    }


