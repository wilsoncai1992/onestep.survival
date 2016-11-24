#' One-step TMLE estimator for survival curve (No censoring)
#'
#' options to ADD:
#' SL.formula: the covariates to include in SL
#'
#' @param dat data.frame with columns T, A, W. All columns with character "W" will be treated as baseline covariates.
#' @param dW binary input vector specifying dynamic treatment (as a function output of W)
#' @param g.lib SuperLearner library for fitting treatment regression
#' @param ht.SL.Lib SuperLearner library for fitting conditional hazard regression
#' @param ... additional options for plotting initial fit curve
#' @param epsilon.step step size for one-step recursion
#' @param max.iter maximal number of recursion for one-step
#' @param tol tolerance for optimization
#' @param verbose to plot the initial fit curve and the objective function value during optimzation
#'
#' @return Psi.hat vector of survival curve under intervention
#' @return T.uniq vector of time points where Psi.hat gets values (have same length as Psi.hat)
#' @return params list of meta-information of estimation
#' @return variables list of data summary
#' @return initial_fit list of initial fit (hazard, g_1, Delta)
#' @export
#'
#' @examples
#' @import dplyr
#' @importFrom plyr rename
#' @import survtmle
#' @import abind
#' @import SuperLearner
surv.one.step.complete <- function(dat,
                                   dW,
                                   g.lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
                                   ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                                   ...,
                                   epsilon.step = 1e-5, # the step size for one-step recursion
                                   max.iter = 1e3, # maximal number of recursion for one-step
                                   tol = 1/nrow(dat), # tolerance for optimization
                                   verbose = TRUE) {
    # ================================================================================================
	# preparation
	# ================================================================================================

	# remove the rows with death time = 0, i.e. who die immediately
	dW <- dW[dat$T.tilde != 0]
	dat <- dat[dat$T.tilde != 0,]

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
	gHatSL <- SuperLearner(Y=dat$A, X=W, SL.library=g.lib, family="binomial")
	# g.hat for each observation
	g.fitted <- gHatSL$SL.predict
	# ================================================================================================
	# conditional hazard (by SL)
	# ================================================================================================
	message('estimating conditional hazard')
	T.uniq <- sort(unique(dat$T.tilde))
	T.max <- max(T.uniq)


	h.hat.t <- haz_SL_wrapper(dat = dat, T.uniq = T.uniq)
	# h.hat at all time t=[0,t.max]
	h.hat.t_full <- as.matrix(h.hat.t$out_haz_full)
	# h.hat at observed unique time t = T.grid
	h.hat.t <- as.matrix(h.hat.t$out_haz)
	# ================================================================================================
	# Qn.A1.t
	# ================================================================================================
	Qn.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

	# ------------------------------------------------------------------------------------
	# compute cumulative hazard
	# ------------------------------------------------------------------------------------
	# cum-product approach (2016-10-05)
	Qn.A1.t_full <- matrix(NA, nrow = n.data, ncol = ncol(h.hat.t_full))
	for (it in 1:n.data) {
	    Qn.A1.t_full[it,] <- cumprod(1 - h.hat.t_full[it,])
	}
	Qn.A1.t <- Qn.A1.t_full[,T.uniq]

	# plot initial fit
	if (verbose) lines(colMeans(Qn.A1.t) ~ T.uniq, ...)
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
	I.A.dW <- dat$A == dW


	D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

	for (it.n in 1:n.data) {
		Y.vec <- create.Y.t.vec(Time = dat$T.tilde[it.n], t.vec = T.uniq)
		temp <- Y.vec - Qn.A1.t[it.n,]
		D1 <- temp / g.fitted[it.n] * I.A.dW[it.n]
		D1.A1 <- temp / g.fitted[it.n] # also update the samples without A = 1
		# D1 matrix
		D1.t[it.n,] <- D1
		D1.A1.t[it.n,] <- D1.A1
	}

	# ================================================================================================
	# Pn.D1: efficient IC average
	# ================================================================================================
	# Pn.D1 vector
	Pn.D1.t <- colMeans(D1.t)
	# ================================================================================================
	# update
	# ================================================================================================
	message('targeting')
	stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq) # 10-17

	update.tensor <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	iter.count <- 0
	stopping.prev <- Inf

	# while ((stopping.criteria >= tol) & (iter.count <= max.iter)) { # ORGINAL
    while ((stopping.criteria >= tol) & (iter.count <= max.iter) & ((stopping.prev - stopping.criteria) >= max(-tol, -1e-5))) { #WILSON: TEMPORARY
		if(verbose) print(stopping.criteria)
		# =============================================================================
		# update the qn
		# ------------------------------------------------------------------------
		# vectorized
		update.mat <- compute.update(D1.t.func.prev = D1.t,
	                                 Pn.D1.func.prev = Pn.D1.t,
	                                 dat = dat,
	                                 T.uniq = T.uniq,
	                                 W_names = W_names,
	                                 dW = dW)
		# ------------------------------------------------------------------------
		update.tensor <- update.tensor + update.mat

		# intergrand <- rowSums(update.tensor)
		# intergrand <- apply(update.tensor, c(1,2), sum)
		intergrand <- update.tensor
		intergrand[is.na(intergrand)] <- 0
		qn.current <- qn.A1.t * exp(epsilon.step * intergrand)

		# normalize the updated qn
		norm.factor <- compute.step.cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf)[,1] #09-06
		qn.current[norm.factor > 1,] <- qn.current[norm.factor > 1,] / norm.factor[norm.factor > 1] #09-06

		# if some qn becomes all zero, prevent NA exisitence
		qn.current[is.na(qn.current)] <- 0
		# =============================================================================
		# compute new Qn

		Qn.current <- compute.step.cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf) # 2016-09-06
		cdf_offset <- 1 - Qn.current[,1] # 2016-09-06
		Qn.current <- Qn.current + cdf_offset # 2016-09-06

		# Qn.current <- apply(qn.current, 1, function(x) compute.step.cdf(pdf.vec = x, t.vec = T.uniq, start = Inf))
		# Qn.current <- t(Qn.current)

		# check error
		# all.equal(compute.step.cdf(pdf.vec = qn.current[1,], t.vec = T.uniq, start = Inf), Qn.current[1,])
		# ------------------------------------------------
		# print('New Qn')

		# < 2016-09-06
		# Qn.current <- matrix(NA, nrow = n.data, ncol = length(T.uniq))
		# for (it.n in 1:n.data) {
			# Qn.current[it.n,] <- rev(cumsum( rev(qn.current[it.n,]) ))
		# }
		# ------------------------------------------------

		# print('New D1')
		# compute new D1
		D1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
		D1.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
		for (it.n in 1:n.data) {
			Y.vec <- create.Y.t.vec(Time = dat$T.tilde[it.n], t.vec = T.uniq)
			temp <- Y.vec - Qn.current[it.n,]
			D1 <- temp / g.fitted[it.n] * I.A.dW[it.n]
			D1.A1 <- temp / g.fitted[it.n] # also update the samples without A = 1
			# D1 matrix
			D1.t[it.n,] <- D1
			D1.A1.t[it.n,] <- D1.A1
		}
		# compute new Pn.D1
		Pn.D1.t <- colMeans(D1.t)
		# ================================================================================================
		# previous stopping criteria
		stopping.prev <- stopping.criteria
		# new stopping criteria
		stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq)
		iter.count <- iter.count + 1

		########################################################################
		# FOR DEBUG ONLY
		# if (TRUE) {
		# 	# ------------------------------------------------------------
		# 	# q <- seq(0,10,.1)
		# 	## truesurvExp <- 1 - pexp(q, rate = 1)
		# 	# truesurvExp <- 1 - pexp(q, rate = .5)
		# 	# plot(round(q*100,0), truesurvExp, type="l", cex=0.2, col = 'red', main = paste('l2 error =', stopping.criteria))
		#
		# 	# library(survival)
		# 	# n.data <- nrow(dat)
		# 	# km.fit <- survfit(Surv(T,rep(1, n.data)) ~ A, data = dat)
		# 	# lines(km.fit)
		# 	# ------------------------------------------------------------
		# 	# Psi.hat <- colMeans(Qn.current[dat$A==dW,])
		#
		# 	# Psi.hat <- colMeans(Qn.current[dat$A==dW & dat$W==0,])
		# 	# ------------------------------------------------------------
		# 	# Q_weighted <- Qn.current/g.fitted[,1] # 09-18: inverse weight by propensity score
		# 	# Q_weighted[dat$A!=dW,] <- 0
		# 	# Psi.hat <- colMeans(Q_weighted) # 09-18: inverse weight by propensity score
		# 	# ------------------------------------------------------------
		#     # 10-06: update all subjects with same W strata
		# 	Psi.hat <- colMeans(Qn.current)
		# 	# ------------------------------------------------------------
		#
		# 	lines(Psi.hat ~ T.uniq, type = 'l', col = 'blue', lwd = .1)
		# 	# ------------------------------------------------------------
		# 	# legend('topright', lty=1, legend = c('true', 'KM', 'one-step'), col=c('red', 'black', 'blue'))
		# }
		########################################################################
		if (iter.count == max.iter) {
			warning('Max Iter count reached, stop iteration.')
		}
	}

	if (!exists('Qn.current')) {
	    # if the iteration immediately converge
	    message('converge suddenly!')
	    Qn.current <- Qn.A1.t
	    Psi.hat <- colMeans(Qn.current)
	}
	# ================================================================================================
	# compute the target parameter
	# ================================================================================================
	# return the mean of those with observed A == dW
	Psi.hat <- colMeans(Qn.current)
	# --------------------------------------------------
	variables <- list(T.uniq = T.uniq)
	params <- list(stopping.criteria = stopping.criteria, epsilon.step = epsilon.step, iter.count = iter.count, max.iter = max.iter, dat = dat)
	initial_fit <- list(h.hat.t = h.hat.t, Qn.A1.t = Qn.A1.t, qn.A1.t = qn.A1.t)
	to.return <- list(Psi.hat = Psi.hat,
	                  T.uniq = T.uniq,
	                  params = params,
	                  variables = variables,
	                  initial_fit = initial_fit)
	class(to.return) <- 'onestep.surv'
	return(to.return)
}
