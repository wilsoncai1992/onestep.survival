compute.step.cdf <- function(pdf.vec, t.vec, start = -Inf) {
		interval.size <- diff(t.vec)
		interval.size <- c(0, interval.size)
		mass.by.interval <- pdf.vec * interval.size
	if(is.infinite(start) & (start < 0)){
		# start from -Inf
		cdf.by.interval <- cumsum(mass.by.interval)
	}else{
		# start from +Inf
		cdf.by.interval <- rev(cumsum(rev(mass.by.interval)))
	}
	
	return(cdf.by.interval)
}


# all.equal(cdf.by.interval, cum.haz.vec)
