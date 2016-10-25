# compute one-step survival curve estimator
# 2016-10-05
# input:
# dat = data.frame with columns T, A, W
# dW = binary input vector as a function output of W
#
# g.lib <- library for g SuperLearner
# ht.SL.Lib <- library for hazard SuperLearner
# verbose <- to plot the initial fit curve and the optimization objective function

# options to ADD:
# SL.formula: the covariates to include in SL


# output:
# Psi.hat <- S_dt function value
# T.uniq <- the jump point of the S_dt
# NOTE: both vector have same length

surv.one.step <- function(dat,
                          dW,
                          col = 'purple', # the col for the verbose curve
                          lty = 1, # col for the verbose curve
                          epsilon.step = 1e-5, # the step size for one-step recursion
                          max.iter = 1e3, # maximal number of recursion for one-step
                          tol = 1/nrow(dat), # tolerance for optimization
                          g.lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
                          ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                          verbose = TRUE
                          ) {
	# ================================================================================================
	# preparation
	# ================================================================================================
	library(dplyr)
	# source('./eval.l2.fast.R')
	source('./eval.l2.step.R')

	# remove the rows with death time = 0, i.e. who die immediately
	# remove the time = T.max, where the censoring probability becomes too large (G(t_|A,W) -> 0) # 10-23
	dW <- dW[dat$T != 0 & dat$T != max(dat$T)]
	dat <- dat[dat$T != 0 & dat$T != max(dat$T),]
    
	
	n.data <- nrow(dat)

	W_names <- grep('W', colnames(dat), value = TRUE)
	W <- dat[,W_names]
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

	gHatSL <- SuperLearner(Y=dat$A, X=W, SL.library=g.lib, family="binomial")
	# g.hat for each observation
	g.fitted <- gHatSL$SL.predict
	# ================================================================================================
	# conditional hazard (by SL)
	# ================================================================================================
	message('estimating conditional hazard')
	T.uniq <- sort(unique(dat$T))
	T.max <- max(T.uniq)

	# source('./ht.SL.R')
	source('./hazard_SL_David/haz_SL_wrapper.R')

	h.hat.t <- haz_SL_wrapper(dat = dat, T.uniq = T.uniq)
	# h.hat at all time t=[0,t.max]
	h.hat.t_full <- as.matrix(h.hat.t$out_haz_full)
	# h.hat at observed unique time t = T.grid
	h.hat.t <- as.matrix(h.hat.t$out_haz)
	# ================================================================================================
	# estimate censoring G(A|W)
	# ================================================================================================
	message('estimating censoring')
	source('./hazard_SL_David/censor_SL_wrapper.R')
	Delta.hat.t <- censor_SL_wrapper(dat = dat, T.uniq = T.uniq,
	                                 Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"))

	# Delta.hat at all time t=[0,t.max]
	Delta.hat.t_full <- as.matrix(Delta.hat.t$out_censor_full)
	# Delta.hat at observed unique time t = T.grid
	Delta.hat.t <- as.matrix(Delta.hat.t$out_censor)
	# ================================================================================================
	# Qn.A1.t
	# ================================================================================================
	source('./gen.step.func.R')
	source('./compute.step.cdf_vec.R')
	Qn.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

	# ------------------------------------------------------------------------------------
	# compute cumulative hazard
	# ------------------------------------------------------------------------------------
	# # integrate approach (2016-09-06)
	# cum.haz.vec <- compute.step.cdf(pdf.mat = h.hat.t, t.vec = T.uniq, start = -Inf)
	# Qn.A1 <- exp( -cum.haz.vec )
	# Qn.A1.t <- Qn.A1

	# cum-product approach (2016-10-05)
	Qn.A1.t_full <- matrix(NA, nrow = n.data, ncol = ncol(h.hat.t_full))
	for (it in 1:n.data) {
	    Qn.A1.t_full[it,] <- cumprod(1 - h.hat.t_full[it,])
	}
	Qn.A1.t <- Qn.A1.t_full[,T.uniq]

	# plot initial fit
	if (verbose) lines(colMeans(Qn.A1.t) ~ T.uniq, col = col, lty = lty)

	# for (it.n in 1:n.data) {
	# 	# print(it.n)
	# 	cum.haz.vec <- rep(NA, length(T.uniq))

	# 	# integrate approach
	#
	# 	# cum.haz.vec <- compute.step.cdf(pdf.vec = h.hat.t[it.n,], t.vec = T.uniq, start = -Inf)
	# 	# ------------------------------------------------------------------------------------
	# 	# cumsum approach
	#
	# 	# cum.haz.vec <- cumsum(h.hat.t[it.n,])
	#
	# 	# haz.func <- gen.step.func(c(h.hat.t[it.n,],0), T.uniq)
	# 	# for (it1 in 1:length(T.uniq)) {
	# 	# 	T.it <- T.uniq[it1]
	# 	# 	# cum.haz <- integrate(haz.func, lower = 0, upper = T.it, subdivisions = 1e5)$value
	# 	# 	cum.haz <- integrate(haz.func, lower = 0, upper = T.it, subdivisions = 1e5, rel.tol = .Machine$double.eps^0.2)$value
	# 	#
	# 	# 	# WILSON: check precision
	# 	# 	cum.haz.vec[it1] <- cum.haz
	# 	# }
	# 	# ------------------------------------------------------------------------------------
	# 	# plot(cum.haz~T.uniq)
	# 	Qn.A1 <- exp( -cum.haz.vec )
	# 	# plot(Qn.A1.t~T.uniq)
	# 	Qn.A1.t[it.n,] <- Qn.A1
	# }

	# plot(Qn.A1.t[20,] ~ T.uniq, type = 'l', xlim=c(0,1e3), ylim = c(0,1))
	# ================================================================================================
	# Gn.A1.t
	# ================================================================================================
	Gn.A1.t <- matrix(NA, nrow = n.data, ncol = length(T.uniq))
	Gn.A1.t_full <- matrix(NA, nrow = n.data, ncol = ncol(Delta.hat.t_full))
	for (it in 1:n.data) {
	    Gn.A1.t_full[it,] <- cumprod(1 - Delta.hat.t_full[it,])
	}
	# G(t_|A,W) = 1 at t = 1, because no censoring at time 0
	Gn.A1.t_full <- cbind(1, Gn.A1.t_full)
	Gn.A1.t_full <- Gn.A1.t_full[, -ncol(Gn.A1.t_full)]
	Gn.A1.t <- Gn.A1.t_full[,T.uniq]

	# Gn.A1.t <- matrix(1, nrow = n.data, ncol = length(T.uniq))
	# Gn.A1.t_full <- matrix(1, nrow = n.data, ncol = T.max)
	
	# plot initial fit
	if (verbose) lines(colMeans(Gn.A1.t) ~ T.uniq, col = 'yellow', lty = lty)
	# ================================================================================================
	# qn.A1.t
	# ================================================================================================
	# WILSON: rewrite in sweep?
	qn.A1.t_full <- matrix(0, nrow = n.data, ncol = ncol(Qn.A1.t_full))
	for (it.n in 1:n.data) {
		qn.A1.t_full[it.n,] <- h.hat.t_full[it.n,] * Qn.A1.t_full[it.n,]
	}
	qn.A1.t <- qn.A1.t_full[,T.uniq]

	# ================================================================================================
	# D1.t: calculate IC
	# D1.A1.t: calculate IC under intervention
	# ================================================================================================
	source('./create.Y.t.vec_censor.R')

	# compute_IC <- function(dat, dW, T.uniq, h.hat.t, g.fitted, Gn.A1.t, Qn.A1.t) {
	#     I.A.dW <- dat$A == dW
	#     n.data <- nrow(dat)
	# 
	#     D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	#     D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	# 
	#     for (it.n in 1:n.data) {
	# 
	#         t_Delta1.vec <- create.Y.t.vec_Delta1(Time = dat$T[it.n], Delta = dat$delta[it.n], t.vec = T.uniq)
	#         t.vec <- create.Y.t.vec(Time = dat$T[it.n], t.vec = T.uniq)
	#         alpha2 <- (t_Delta1.vec - t.vec * h.hat.t[it.n,])
	# 
	#         alpha1 <- -I.A.dW[it.n]/g.fitted[it.n]/Gn.A1.t[it.n,]/Qn.A1.t[it.n,]
	#         alpha1_A1 <- -1/g.fitted[it.n]/Gn.A1.t[it.n,]/Qn.A1.t[it.n,]
	#         
	#         not_complete <- alpha1 * alpha2
	#         not_complete_A1 <- alpha1_A1 * alpha2
	#         # D1 matrix
	#         D1.t[it.n, ] <- cumsum(not_complete) * Qn.A1.t[it.n,] # complete influence curve
	#         D1.A1.t[it.n, ] <- cumsum(not_complete_A1) * Qn.A1.t[it.n,] # also update those A = 0.
	#     }
	#     
	#     # turn unstable results to 0
	#     D1.t[is.na(D1.t)] <- 0
	#     D1.A1.t[is.na(D1.A1.t)] <- 0
	#     
	#     return(list(D1.t = D1.t,
	#                 D1.A1.t = D1.A1.t))
	# }

	
	compute_IC <- function(dat, dW, T.uniq, h.hat.t_full, g.fitted, Gn.A1.t_full, Qn.A1.t, Qn.A1.t_full) {
	    I.A.dW <- dat$A == dW
	    n.data <- nrow(dat)
	    
	    D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	    D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	    
	    for (it.n in 1:n.data) {
	        
	        t_Delta1.vec <- create.Y.t.vec_Delta1(Time = dat$T[it.n], Delta = dat$delta[it.n], t.vec = 1:max(T.uniq))
	        t.vec <- create.Y.t.vec(Time = dat$T[it.n], t.vec = 1:max(T.uniq))
	        alpha2 <- (t_Delta1.vec - t.vec * h.hat.t_full[it.n,])
	        
	        alpha1 <- -I.A.dW[it.n]/g.fitted[it.n]/Gn.A1.t_full[it.n,]/Qn.A1.t_full[it.n,]
	        alpha1_A1 <- -1/g.fitted[it.n]/Gn.A1.t_full[it.n,]/Qn.A1.t_full[it.n,]
	        
	        not_complete <- alpha1 * alpha2
	        not_complete_A1 <- alpha1_A1 * alpha2
	        # D1 matrix
	        D1.t[it.n, ] <- cumsum(not_complete)[T.uniq] * Qn.A1.t[it.n,] # complete influence curve
	        D1.A1.t[it.n, ] <- cumsum(not_complete_A1)[T.uniq] * Qn.A1.t[it.n,] # also update those A = 0.
	    }
	    
	    # turn unstable results to 0
	    D1.t[is.na(D1.t)] <- 0
	    D1.A1.t[is.na(D1.A1.t)] <- 0
	    
	    return(list(D1.t = D1.t,
	                D1.A1.t = D1.A1.t))
	}
	
	initial_IC <- compute_IC(dat = dat,
	                         dW = dW,
	                         T.uniq = T.uniq,
	                         h.hat.t_full = h.hat.t_full,
	                         g.fitted = g.fitted,
	                         Gn.A1.t_full = Gn.A1.t_full,
	                         Qn.A1.t = Qn.A1.t,
	                         Qn.A1.t_full = Qn.A1.t_full)
	D1.t <- initial_IC$D1.t
	D1.A1.t <- initial_IC$D1.A1.t
	# ================================================================================================
	# Pn.D1: efficient IC average
	# ================================================================================================
	# Pn.D1 vector
	Pn.D1.t <- colMeans(D1.t)
	# Pn.D1.t.func <- gen.step.func(y.vec = c(0, Pn.D1.t), t.vec = T.uniq)

	# D1.t.func <- gen.step.func(y.vec = c(1, D1.t[it.n,]), t.vec = T.uniq)
	## curve(Qn.A1.t.func, from = 0, to = 400)
	## curve(Pn.D1.t.func, from = 0, to = 400)
	# ================================================================================================
	# update
	# ================================================================================================
	source('./compute_update.R')
	# ================================================================================================
	# stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))
	stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq) # 10-17

    # update.tensor <- array(0, dim = c(n.data ,length(T.uniq) ,2))
	update.tensor <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	iter.count <- 0
	stopping.prev <- Inf

	library(abind)
	# while ((stopping.criteria >= tol) & (iter.count <= max.iter)) { # ORGINAL
    while ((stopping.criteria >= tol) & (iter.count <= max.iter) & ((stopping.prev - stopping.criteria) >= max(-tol, -1e-5))) { #WILSON: TEMPORARY
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
		update.mat <- compute.update(D1.t.func.prev = D1.t,
	                                 Pn.D1.func.prev = Pn.D1.t,
	                                 dat = dat,
	                                 T.uniq = T.uniq,
	                                 W_names = W_names,
	                                 dW = dW)
		# ------------------------------------------------------------------------
		# update.tensor <- abind(update.tensor, update.mat, rev.along = 1)
		update.tensor <- update.tensor + update.mat

		# intergrand <- rowSums(update.tensor)
		# intergrand <- apply(update.tensor, c(1,2), sum)
		intergrand <- update.tensor
		intergrand[is.na(intergrand)] <- 0
		qn.current <- qn.A1.t * exp(epsilon.step * intergrand)
		qn.current_full <- qn.A1.t_full * exp(epsilon.step * replicate(T.max, intergrand[,1])) #10-23

		# normalize the updated qn
		# norm.factor <- compute.step.cdf(qn.current, t.vec = T.uniq, start = -Inf)
		# norm.factor <- norm.factor[,ncol(norm.factor)]

		norm.factor <- compute.step.cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf)[,1] #09-06
		qn.current[norm.factor > 1,] <- qn.current[norm.factor > 1,] / norm.factor[norm.factor > 1] #09-06
		qn.current_full[norm.factor > 1,] <- qn.current_full[norm.factor > 1,] / norm.factor[norm.factor > 1] #10-23
		
		# norm.factor <- rowSums(qn.current)
		# qn.current <- qn.current / norm.factor
		# if some qn becomes all zero, prevent NA exisitence
		qn.current[is.na(qn.current)] <- 0
		qn.current_full[is.na(qn.current_full)] <- 0 #10-23
		# =============================================================================
		# compute new Qn

		Qn.current <- compute.step.cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf) # 2016-09-06
		cdf_offset <- 1 - Qn.current[,1] # 2016-09-06
		Qn.current <- Qn.current + cdf_offset # 2016-09-06

		Qn.current_full <- compute.step.cdf(pdf.mat = qn.current_full, t.vec = 1:max(T.uniq), start = Inf) # 10-23
		cdf_offset <- 1 - Qn.current_full[,1] # 10-23
		Qn.current_full <- Qn.current_full + cdf_offset # 10-23
		
		# check error
		# all.equal(compute.step.cdf(pdf.vec = qn.current[1,], t.vec = T.uniq, start = Inf), Qn.current[1,])
		# ------------------------------------------------
		# print('New Qn')

		# < 2016-09-06
		# Qn.current <- matrix(NA, nrow = n.data, ncol = length(T.uniq))
		# for (it.n in 1:n.data) {
			# Qn.current[it.n,] <- rev(cumsum( rev(qn.current[it.n,]) ))
		# }
		# =============================================================================
		# compute new h_t
		h.hat.t_full_current <- matrix(0, nrow = n.data, ncol = max(T.uniq))
		for (it.n in 1:n.data) {
		    h.hat.t_full_current[it.n, ] <- qn.current_full[it.n, ] / Qn.current_full[it.n,]
		}
		# ------------------------------------------------
		# print('New D1')
		# compute new D1

		# updated_IC <- compute_IC(dat = dat,
		#                          dW = dW,
		#                          T.uniq = T.uniq,
		#                          h.hat.t = h.t_current,
		#                          g.fitted = g.fitted,
		#                          Gn.A1.t = Gn.A1.t,
		#                          Qn.A1.t = Qn.current)

		updated_IC <- compute_IC(dat = dat,
		                         dW = dW,
		                         T.uniq = T.uniq,
		                         h.hat.t_full = h.hat.t_full_current,
		                         g.fitted = g.fitted,
		                         Gn.A1.t_full = Gn.A1.t_full,
		                         Qn.A1.t = Qn.current,
		                         Qn.A1.t_full = Qn.current_full)
		
		# compute new D1
		# D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
		# D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
		# for (it.n in 1:n.data) {
		# 	Y.vec <- create.Y.t.vec(Time = dat$T[it.n], t.vec = T.uniq)
		# 	temp <- Y.vec - Qn.current[it.n,]
		# 	D1 <- temp / g.fitted[it.n] * I.A.dW[it.n]
		# 	D1.A1 <- temp / g.fitted[it.n] # also update the samples without A = 1
		# 	# D1 matrix
		# 	D1.t[it.n,] <- D1
		# 	D1.A1.t[it.n,] <- D1.A1
		# }
		D1.t <- updated_IC$D1.t
		D1.A1.t <- updated_IC$D1.A1.t
		
		# compute new Pn.D1
		Pn.D1.t <- colMeans(D1.t)


		# Pn.D1.t.func <- gen.step.func(y.vec = c(0, Pn.D1.t), t.vec = T.uniq)
		#
		# D1.t.func <- gen.step.func(y.vec = c(1, D1.t[it.n,]), t.vec = T.uniq)
		# ================================================================================================
		# previous stopping criteria
		stopping.prev <- stopping.criteria
		# new stopping criteria
		# stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))
		stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq)
		# interation count += 1
		iter.count <- iter.count + 1

		########################################################################
		# FOR DEBUG ONLY
		if (TRUE) {
			# ------------------------------------------------------------
			# q <- seq(0,10,.1)
			## truesurvExp <- 1 - pexp(q, rate = 1)
			# truesurvExp <- 1 - pexp(q, rate = .5)
			# plot(round(q*100,0), truesurvExp, type="l", cex=0.2, col = 'red', main = paste('l2 error =', stopping.criteria))

			# library(survival)
			# n.data <- nrow(dat)
			# km.fit <- survfit(Surv(T,rep(1, n.data)) ~ A, data = dat)
			# lines(km.fit)
			# ------------------------------------------------------------
			# Psi.hat <- colMeans(Qn.current[dat$A==dW,])

			# Psi.hat <- colMeans(Qn.current[dat$A==dW & dat$W==0,])
			# ------------------------------------------------------------
			# Q_weighted <- Qn.current/g.fitted[,1] # 09-18: inverse weight by propensity score
			# Q_weighted[dat$A!=dW,] <- 0
			# Psi.hat <- colMeans(Q_weighted) # 09-18: inverse weight by propensity score
			# ------------------------------------------------------------
		    # 10-06: update all subjects with same W strata
			Psi.hat <- colMeans(Qn.current)
			# ------------------------------------------------------------

			lines(Psi.hat ~ T.uniq, type = 'l', col = 'blue', lwd = .1)
			# ------------------------------------------------------------
			# legend('topright', lty=1, legend = c('true', 'KM', 'one-step'), col=c('red', 'black', 'blue'))
		}
		########################################################################
		if (iter.count == max.iter) {
			warning('Max Iter count reached, stop iteration.')
		}
	}

	# ================================================================================================
	# compute the target parameter
	# ================================================================================================
	# return the mean of those with observed A == dW
	# Qn.current <- Qn.current[dat$A==dW,]
	# Psi.hat <- colMeans(Qn.current)
	# --------------------------------------------------
	# 09-18: inverse weight by propensity score
	# Q_weighted <- Qn.current/g.fitted[,1]
	# Q_weighted[dat$A!=dW,] <- 0
	# Psi.hat <- colMeans(Q_weighted)

	# --------------------------------------------------
	# 10-06: update all subjects
	Psi.hat <- colMeans(Qn.current)
	# --------------------------------------------------
	variables <- list(T.uniq)
	params <- list(stopping.criteria, epsilon.step, iter.count, max.iter)
	initial_fit <- list(h.hat.t, Qn.A1.t, qn.A1.t)
	to.return <- list(Psi.hat = Psi.hat,
	                  T.uniq = T.uniq,
	                  params = params,
	                  variables = variables,
	                  initial_fit = initial_fit)
	return(to.return)
}