# OLD: 2016-05-30
surv.one.step <- function(dat, dW) {
	# ================================================================================================
	# preparation
	# ================================================================================================
	# source('./eval.l2.fast.R')
	source('./eval.l2.step.R')
	epsilon.step <- 1e-2
	# epsilon.step <- 1e-2
	max.iter <- 1e2
	n.data <- nrow(dat)
	tol <- 1/n.data
	
	W <- dat$W
	W <- as.data.frame(W)
	# ================================================================================================
	# input validation
	# ================================================================================================
	if (length(dW) != n.data) {
		stop('The length of input dW is not same as the sample size!')
	}
	# ================================================================================================
	# estimate g(A|W)
	# ================================================================================================
	library(SuperLearner)
	SL.library<- c("SL.glm", "SL.step", "SL.glm.interaction" )
	gHatSL <- SuperLearner(Y=dat$A, X=W, SL.library=SL.library, family="binomial")
	# g.hat for each observation
	g.fitted <- gHatSL$SL.predict
	# ================================================================================================
	# conditional hazard (by SL)
	# ================================================================================================
	print('estimating conditional hazard')
	T.uniq <- sort(unique(dat$T))
	h.hat.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

	source('./ht.SL.R')
	for (it1 in 1:length(T.uniq)) {
		print(it1)
		T.it <- T.uniq[it1]
		I.t <- (dat$T == T.it) + 0
		# ----------------------------------------------------------------------------------------
		# SL for h.n
		# ----------------------------------------------------------------------------------------
		hn.A1 <- ht.SL(dat, T.it, I.t)
		# ----------------------------------------------------------------------------------------
		# X.mat <- dat[,c('A', 'W')]
		# predict counterfactual
		# X.mat.new <- X.mat
		# X.mat.new[,1] <- 1
		# hHat.SL <- SuperLearner(Y=I.t, X=X.mat, newX = X.mat.new,
														# SL.library=SL.library, family="binomial")
		# hn.A1 <- hHat.SL$SL.predict

		# ----------------------------------------------------------------------------------------
		h.hat.t[,it1] <- hn.A1
	}
	# plot(h.hat.t[1,] ~ T.uniq, xlim = c(0,1e3), ylim = c(0,10))
	# ================================================================================================
	# TEMP: conditional hazard by nelson aalen
	# ================================================================================================
	# source('./na.haz.R')
	# play <- na.haz(dat)
	# 
	# h.hat.t <- rep(0, length(T.uniq))
	# ind <- match(play[[2]], T.uniq)
	# h.hat.t[ind] <- play[[1]]
	# rep.row<-function(x,n){
	# 	matrix(rep(x,each=n),nrow=n)
	# }
	# h.hat.t <- rep.row(h.hat.t, n.data)
	# ================================================================================================
	# Qn.A1.t
	# ================================================================================================
	source('./gen.step.func.R')
	# source('./compute.step.cdf.R')
	source('./compute.step.cdf_vec.R')
	Qn.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	for (it.n in 1:n.data) {
		# print(it.n)
		cum.haz.vec <- rep(NA, length(T.uniq))
		# ------------------------------------------------------------------------------------
		# compute cumulative hazard
		# ------------------------------------------------------------------------------------
		# integrate approach
		# cum.haz.vec <- compute.step.cdf(pdf.vec = h.hat.t[it.n,], t.vec = T.uniq, start = -Inf)
		# ------------------------------------------------------------------------------------
		# cumsum approach
		cum.haz.vec <- cumsum(h.hat.t[it.n,])
		
		# haz.func <- gen.step.func(c(h.hat.t[it.n,],0), T.uniq)
		# for (it1 in 1:length(T.uniq)) {
		# 	T.it <- T.uniq[it1]
		# 	# cum.haz <- integrate(haz.func, lower = 0, upper = T.it, subdivisions = 1e5)$value
		# 	cum.haz <- integrate(haz.func, lower = 0, upper = T.it, subdivisions = 1e5, rel.tol = .Machine$double.eps^0.2)$value
		#
		# 	# WILSON: check precision
		# 	cum.haz.vec[it1] <- cum.haz
		# }
		# ------------------------------------------------------------------------------------
		# plot(cum.haz~T.uniq)
		Qn.A1 <- exp( -cum.haz.vec )
		# plot(Qn.A1.t~T.uniq)
		Qn.A1.t[it.n,] <- Qn.A1
	}
	# plot(Qn.A1.t[20,] ~ T.uniq)
	# ================================================================================================
	# qn.A1.t
	# ================================================================================================
	qn.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	for (it.n in 1:n.data) {
		qn.A1 <- h.hat.t[it.n,] * Qn.A1.t[it.n,]
		# plot(qn.A1~T.uniq, type = 'l')
		qn.A1.t[it.n,] <- qn.A1
	}
	# ================================================================================================
	# D1.t: calculate IC
	# ================================================================================================
	I.A.dW <- dat$A == dW
	
	source('./create.Y.t.vec.R')
	
	D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	for (it.n in 1:n.data) {
		Y.vec <- create.Y.t.vec(Time = dat$T[it.n], t.vec = T.uniq)
		temp <- Y.vec - Qn.A1.t[it.n,]
		D1 <- temp / g.fitted[it.n] * I.A.dW[it.n]
		# D1 matrix
		D1.t[it.n,] <- D1
	}
	# ================================================================================================
	# Pn.D1: efficient IC average
	# ================================================================================================
	# Pn.D1 vector
	Pn.D1.t <- colMeans(D1.t)
	# Pn.D1.t.func <- gen.step.func(y.vec = c(0, Pn.D1.t), t.vec = T.uniq)
	
	# D1.t.func <- gen.step.func(y.vec = c(1, D1.t[it.n,]), t.vec = T.uniq)
	## curve(Qn.A1.t.func, from = 0, to = 10000)
	# ================================================================================================
	# update
	# ================================================================================================
	compute.update <- function(D1.t.func.prev, Pn.D1.func.prev) {
		# formula on p.30
		# result <- l2.inner.step(Pn.D1.t, D1.t, T.uniq) /
					# sqrt(l2.inner.step(D1.t, D1.t, T.uniq))
		# formula on p.28
		# result <- l2.inner.step(Pn.D1.func.prev, D1.t.func.prev, T.uniq) /
			# sqrt(l2.inner.step(Pn.D1.func.prev, Pn.D1.func.prev, T.uniq))
		
		# WILSON MADE: MAY BE WRONG
		result <- l2.inner.step(abs(Pn.D1.func.prev), D1.t.func.prev, T.uniq) /
			sqrt(l2.inner.step(Pn.D1.func.prev, Pn.D1.func.prev, T.uniq))
		
		# if (is.na(result)) {
		# result <- NA
		# }
		# result <- as.vector(result)
		return(result)
	}
	# ================================================================================================
	T.max <- max(T.uniq)
	
	stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))
	update.seq <- matrix(0, nrow = n.data, ncol = 1)
	iter.count <- 0
	stopping.prev <- Inf
	
	while ((stopping.criteria >= tol) & (iter.count <= max.iter) & (stopping.criteria <= stopping.prev)) {
	# while ((stopping.criteria >= tol) & (iter.count <= max.iter)) { #WILSON: TEMPORARY
		print(stopping.criteria)
		# print('compute update value')
		# =============================================================================
		# update the qn
		# ------------------------------------------------------------------------
		# WILSON: not vectorized
		# update.vec <- rep(0, n.data)
		# for (it.n in 1:n.data) {
		# 	print(it.n)
		# 	## D1.t.func <- gen.step.func(y.vec = c(1, D1.t[it.n,]), t.vec = T.uniq)
		# 	update.val <- compute.update(D1.t.func.prev = D1.t,
		# 															 Pn.D1.func.prev = Pn.D1.t)
		# 	update.vec[it.n] <- update.val
		# }
		# ------------------------------------------------------------------------
		# vectorized
		update.vec <- compute.update(D1.t.func.prev = D1.t,
																 Pn.D1.func.prev = Pn.D1.t)
		# ------------------------------------------------------------------------
		update.seq <- cbind(update.seq, update.vec)
		
		intergrand <- rowSums(update.seq)
		intergrand[is.na(intergrand)] <- 0
		qn.current <- qn.A1.t * exp(epsilon.step * intergrand)
		
		# =============================================================================
		# compute new Qn
		# Qn.current <- apply(qn.current, 1, function(x) compute.step.cdf(pdf.vec = x, t.vec = T.uniq, start = Inf))
		# Qn.current <- t(Qn.current)
		
		# check error
		# all.equal(compute.step.cdf(pdf.vec = qn.current[1,], t.vec = T.uniq, start = Inf), Qn.current[1,])
		# ------------------------------------------------
		# print('New Qn')
		Qn.current <- matrix(NA, nrow = n.data, ncol = length(T.uniq))
		for (it.n in 1:n.data) {
			Qn.current[it.n,] <- rev(cumsum( rev(qn.current[it.n,]) ))
		}
		# ------------------------------------------------
		
		# print('New D1')
		# compute new D1
		D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
		for (it.n in 1:n.data) {
			Y.vec <- create.Y.t.vec(Time = dat$T[it.n], t.vec = T.uniq)
			temp <- Y.vec - Qn.current[it.n,]
			D1 <- temp / g.fitted[it.n] * I.A.dW[it.n]
			# D1 matrix
			D1.t[it.n,] <- D1
		}
		# compute new Pn.D1
		Pn.D1.t <- colMeans(D1.t)
		# Pn.D1.t.func <- gen.step.func(y.vec = c(0, Pn.D1.t), t.vec = T.uniq)
		# 
		# D1.t.func <- gen.step.func(y.vec = c(1, D1.t[it.n,]), t.vec = T.uniq)
		# ================================================================================================
		# previous stopping criteria
		stopping.prev <- stopping.criteria
		# new stopping criteria
		stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))
		# interation count += 1
		iter.count <- iter.count + 1
		
		########################################################################
		# FOR DEBUG ONLY
		if (TRUE) {
			Psi.hat <- colMeans(Qn.current)
			lines(Psi.hat ~ T.uniq, type = 'l', col = 'blue')
		}
		########################################################################
	}
	
	# ================================================================================================
	# compute the target parameter
	# ================================================================================================
	########################################################################
	# FOR DEBUG ONLY
	########################################################################
	# Qn.current[dat$A==0,] <- matrix(0, ncol = length(T.uniq), nrow = sum(dat$A==0))
	# Qn.current <- Qn.current[dat$A==1,]
	########################################################################
	Psi.hat <- colMeans(Qn.current)
	
	to.return <- list(Psi.hat, T.uniq)
	return(to.return)
}
# ================================================================================================
#
# ================================================================================================
# for (play in 1:n.data) {
# 	play <- play + 1
# 	plot(h.hat.t[play,] ~ T.uniq, type = 'l')
# }
# 
# lines(Psi.hat ~ T.uniq, lty=2)
