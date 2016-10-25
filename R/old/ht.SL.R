# a poor-performing hazard estimator 
# 2016-06-01

ht.SL <- function(dat, T.it, I.t) {
	X.mat <- dat[,c('A', 'W')]
	# predict counterfactual
	X.mat.new <- X.mat
	X.mat.new[,1] <- 1
	# estimate at risk probability
	X.mat.2 <- X.mat
	####################################
	# Temporaray!
	# subset for conditional hazard
	is.at.risk <- (dat$T >= T.it)
	X.mat <- X.mat[is.at.risk,]
	I.t <- I.t[is.at.risk]
	####################################
	
	
	# ----------------------------------------------------------------------------------------
	# use glm for faster
	X.mat <- cbind(I.t, X.mat)
	hhat.glm <- glm(I.t ~ A + W, family = 'binomial', data = X.mat)
	hn.A1.numerator <- predict(hhat.glm, newdata = X.mat.new, type = 'response')
	# ----------------------------------------------------------------------------------------	
	# fit denominator: at risk | A,W
	X.mat.2 <- cbind(is.at.risk, X.mat.2, row.names = NULL)
	hhat.glm.2 <- glm(is.at.risk ~ A + W, family = 'binomial', data = X.mat.2)
	hn.A1.denominator <- predict(hhat.glm.2, newdata = X.mat.new, type = 'response')
	# ----------------------------------------------------------------------------------------	
	# calculate ht
	hn.A1 <- hn.A1.numerator / hn.A1.denominator
	return(hn.A1)
}