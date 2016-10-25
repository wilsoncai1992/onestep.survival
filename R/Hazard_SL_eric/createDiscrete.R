# createDiscrete takes a survival time data set and returns an expanded binary process data set

## previous verious:
createDiscrete <- function(time, event, dataX, delta.upper = NULL, n.delta=100) {
	# time is the follow up time
	# event, 1 is event, 0 is censored
	# dataX contains all X variables
	# n.delta is the number of time windows
	n <- length(time)
	if(is.null(delta.upper)) {
		n.delta <- min(n.delta, length(unique(time[event==1])))
		probs.delta <- seq(from=0, to=1, length.out=n.delta)
		delta.upper <- quantile(time[event==1], probs=probs.delta, names=FALSE)
	}
	delta.lower <- c(0, delta.upper[-length(delta.upper)])
	n.delta <- length(delta.upper)
	ID <- seq(n) 
	DATA <- cbind(ID, time, event, dataX)
	long.DATA <- apply(DATA, 2, function(x) rep(x, times=rep(n.delta, times=n)))
	N.delta <- rep(NA, nrow(long.DATA))
	long.DATA <- cbind(long.DATA, delta.lower, delta.upper, N.delta)
	long.DATA <- as.data.frame(long.DATA)
	long.DATA$N.delta <- ifelse(long.DATA$time > long.DATA$delta.upper, 0, ifelse(long.DATA$event==1, ifelse(long.DATA$time <= long.DATA$delta.lower, NA, 1), NA))
	long.DATA <- long.DATA[!is.na(long.DATA$N.delta), ]
	return(long.DATA)	
}


# updated version:
# createDiscrete <- function(time, event, dataX, lengthDelta = 100, deltaUpper = NULL) {
# 	# time is the follow up time
# 	# event, 1 is event, 0 is censored
# 	# dataX contains all X variables as matrix
# 	# lengthDelta is the number of time windows
# 	# alternative to lengthDelta is to give exact cut points by a vector specifying the upper limit, deltaUpper
# 	n <- length(time)
# 	if(is.null(deltaUpper)) {
# 		# number of cut points should be larger than number of unique event times
# 		lengthDelta <- min(lengthDelta, length(unique(time[event==1])))
# 		# use quantiles of unique event times to create cut points
# 		probs.delta <- seq(from=0, to=1, length.out=lengthDelta)
# 		deltaUpper <- quantile(time[event==1], probs=probs.delta, names=FALSE)
# 	}
# 	deltaLower <- c(0, deltaUpper[-length(deltaUpper)])
# 	lengthDelta <- length(deltaUpper)
# 	ID <- seq(n) 
# 	DATA <- cbind(ID, time, event, dataX)
# 	# replicate DATA lengthDelta times
# 	longDATA <- apply(DATA, 2, function(x) rep(x, times=rep.int(lengthDelta, times=n)))
# 	N.delta <- rep.int(NA, times = nrow(longDATA))
# 	longDATA <- cbind(longDATA, deltaLower, deltaUpper, N.delta)
# 	# update N.delta with the appropriate value
# 	# N(delta) = 0 if still at risk
# 	# N(delta) = 1 if have event in [deltaLower, deltaUpper)
# 	# N(delta) = NA if censored or had event
# 	longDATA[, "N.delta"] <- ifelse(longDATA[, "time"] > longDATA[, "deltaUpper"], 0, ifelse(longDATA[, "event"]==1, ifelse(longDATA[, "time"] <= longDATA[, "deltaLower"], NA, 1), NA))
# 	# reduce to observed outcomes
# 	longDATA <- longDATA[!is.na(longDATA[, "N.delta"]), ]
# 	# convert to data.frame?
# 	longDATA <- as.data.frame(longDATA)
# 	return(longDATA)	
# }
