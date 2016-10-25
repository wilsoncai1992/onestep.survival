na.haz <- function(dat) {
	library(mice)
	status <- rep(1, sum(dat$A==1))
	time <- dat$T[dat$A==1]
	time <- sort(time)
	eng <- data.frame(time, status)
	ch <- nelsonaalen(eng, time, status)
	# plot(x = time, y = ch, ylab='Cumulative hazard', xlab='Time')
	
	
	ch <- ch[as.numeric(rownames(unique(data.frame(time)[1])))]
	ha <- diff(c(0,ch))
	# plot(ha ~ unique(time))
	# summary(ha)	
	
	
	return(list(ha, unique(time)))
}
