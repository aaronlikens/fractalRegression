#' Iterated Amplitude Adjusted Fourier Transform
#' 
#' @param signal is a real valued time serires 
#' @param N is the number of desired surrogates. Default is 1
#' @export
iaafft <- function(signal,N=1){
# this function generates surrogates using the iterated amplitude
# adjusted fourier transform discussed in Ihlen & Vereijken, 2010
# and Schreiber & Schmitz, 1996 and gazillion other papers.


	mx = 1000
	x = signal
	ln = length(x)
	amp = Mod(fft(x))
	sgates = matrix(rep(0,ln*N),nrow = ln)
	for (n in 1:N){
		s = sample(ln)
		sgates[,n]<-x[s]
	}
	
	tmp <- sort(x,index.return = TRUE)
	x <- tmp$x
	ind <- tmp$ix
	
	for ( n in 1:N){
		phase_x <- Arg(fft(sgates[,n]))
		nn =1
		conv = 0
		ind_prev=ind
		while ( nn <= mx && conv == 0 ){
			sgates[,n] <- amp*exp(phase_x*1i)
			sgates[,n] <- Re(fft(sgates[,n],inverse = TRUE))
			sgates_sort_ind <- sort(sgates[,n],index.return = TRUE)
			sgates_sort_new <- sort(sgates_sort_ind$ix,index.return=TRUE)
			sgates[,n] <- x[sgates_sort_new$ix]
			ind_new <- sgates_sort_new$ix
			if (all(ind_new == ind_prev)){
				conv <-1
			}else{
				ind_prev <- ind_new
				nn = nn+1;
			}
			phase_x <- Arg(fft(sgates[,n]))
		}
	}
	sgates <- data.frame(Re(sgates))
	return(sgates)
}

 