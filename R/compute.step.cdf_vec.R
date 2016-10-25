# compute cdf of step function
# pdf.mat: if input vector = compute cdf for a step-function pdf
# if input matrix = compute cdf for several step-function pdf with same jump points

# t.vec: unique jump points of step function
# start: -Inf = from left to right
# Inf = from right to left

# 2016-06-01
# VECTORIZED VERSION
compute.step.cdf <- function(pdf.mat, t.vec, start = -Inf) {
		interval.size <- diff(t.vec)
		# interval.size <- c(0, interval.size)
		interval.size <- c(interval.size, 0) # 09-07
		
		# compute the mass
		if(is.matrix(pdf.mat)){
			# if input with multi-sample
			mass.by.interval <- sweep(pdf.mat,MARGIN=2, interval.size, `*`) 
			# multiplies the interval length to each row of the y-values
			# the result is a matrix, each row is a single pdf, and entries are the mass
			
		}else{
			# if input with one-sample
			mass.by.interval <- pdf.mat * interval.size
		}
		
		
	if(is.infinite(start) & (start < 0)){
		# ======================================================================
		# start from -Inf
		if(is.matrix(pdf.mat)){
			# if input with multi-sample
			cdf.by.interval <- t(apply(mass.by.interval, 1, cumsum)) # cumsum of mass for each row, from left to right
		}else{
			# if input with one-sample
			cdf.by.interval <- cumsum(mass.by.interval)
		}
	}else{
		# ======================================================================
		# start from +Inf
		if(is.matrix(pdf.mat)){
			# if input with multi-sample
			cdf.by.interval <- t(apply(mass.by.interval, 1, function(obj) rev(cumsum(rev(obj))) ) )
		}else{
			# if input with one-sample
			cdf.by.interval <- rev(cumsum(rev(mass.by.interval)))
		}
	}
		# ======================================================================
	return(cdf.by.interval)
}


# all.equal(cdf.by.interval, cum.haz.vec)
