library(mice)
require(MASS)

leuk$status <- 1  ## no censoring occurs in leuk data (MASS)
ch <- nelsonaalen(leuk, time, status)
plot(x = leuk$time, y = ch, ylab='Cumulative hazard', xlab='Time')

### See example on <a href="http://www.engineeredsoftware.com/lmar/pe_cum_hazard_function.htm<br />
# time" title="http://www.engineeredsoftware.com/lmar/pe_cum_hazard_function.htm<br />
	# time">http://www.engineeredsoftware.com/lmar/pe_cum_hazard_function.htm<br />
time <- c(43, 67, 92, 94, 149, rep(149,7))
status <- c(rep(1,5),rep(0,7))
eng <- data.frame(time, status)
ch <- nelsonaalen(eng, time, status)
plot(x = time, y = ch, ylab='Cumulative hazard', xlab='Time')

diff(ch)
time

# ==========
# status <- rep(1, nrow(dat))
# time <- dat$T

na.haz <- function(dat) {
	library(mice)
	status <- rep(1, sum(dat$A==1))
	time <- dat$T[dat$A==1]
	time <- sort(time)
	eng <- data.frame(time, status)
	ch <- nelsonaalen(eng, time, status)
	plot(x = time, y = ch, ylab='Cumulative hazard', xlab='Time')
	
	
	ch <- ch[as.numeric(rownames(unique(data.frame(time)[1])))]
	ha <- diff(c(0,ch))
	# plot(ha ~ unique(time))
	# summary(ha)	
	return(ha)
}

# =================================================================================
fit1 <- muhaz(time, status)
plot(fit1)
summary(fit1)
# =================================================================================
# Nelson aalen estimator
# =================================================================================
calcna = function(time, event) {
	na.fit = survfit(coxph(Surv(time,event)~1), type="aalen")
	jumps = c(0, na.fit$time, max(time))
	# need to be careful at the beginning and end
	surv = c(1, na.fit$surv, na.fit$surv[length(na.fit$surv)])
	
	# apply appropriate transformation
	neglogsurv = -log(surv)   
	
	# create placeholder of correct length
	naest = numeric(length(time))  
	for (i in 2:length(jumps)) {
		naest[which(time>=jumps[i-1] & time<=jumps[i])] = 
			neglogsurv[i-1]   # snag the appropriate value
	}
	return(naest)
}

newna = calcna(time, status)
yo <- cbind(time, newna)
plot(newna ~ time, data = yo)
