#' compute one-step survival curve estimator
#'
#' options to ADD:
#' SL.formula: the covariates to include in SL
#'
#' @param dat data.frame with columns T, A, C, W. All columns with character "W" will be treated as baseline covariates.
#' @param dW binary input vector specifying dynamic treatment (as a function output of W)
#' @param g.lib SuperLearner library for fitting treatment regression
#' @param Delta.SL.Lib SuperLearner library for fitting censoring regression
#' @param ht.SL.Lib SuperLearner library for fitting conditional hazard regression
#' @param ... additional options for plotting initial fit curve
#' @param epsilon.step step size for one-step recursion
#' @param max.iter maximal number of recursion for one-step
#' @param tol  tolerance for optimization
#' @param verbose to plot the initial fit curve and the objective function value during optimzation
#'
#' @return Psi.hat vector of survival curve under intervention
#' @return T.uniq vector of time points where Psi.hat gets values (have same length as Psi.hat)
#' @return params list of meta-information of estimation
#' @return variables list of data summary
#' @return initial_fit list of initial fit (hazard, g_1, Delta)
#'
#' @export
#'
#' @examples
#' # TODO
#' @import dplyr
#' @importFrom plyr rename
#' @import survtmle
#' @import abind
#' @import SuperLearner
surv.one.step <- function(dat,
                          dW,
                          g.lib = c("SL.glm", "SL.step", "SL.glm.interaction"),
                          Delta.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                          ht.SL.Lib = c("SL.mean","SL.glm", "SL.gam", "SL.earth"),
                          ...,
                          epsilon.step = 1e-5,
                          max.iter = 1e3,
                          tol = 1/nrow(dat),
                          verbose = TRUE) {
	# ================================================================================================
	# preparation
	# ================================================================================================
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
	gHatSL <- SuperLearner(Y=dat$A, X=W, SL.library=g.lib, family="binomial")
	# g.hat for each observation
	g.fitted <- gHatSL$SL.predict
	# ================================================================================================
	# conditional hazard (by SL)
	# ================================================================================================
	message('estimating conditional hazard')
	T.uniq <- sort(unique(dat$T))
	T.max <- max(T.uniq)


	h.hat.t <- haz_SL_wrapper(dat = dat, T.uniq = T.uniq)
	# h.hat at all time t=[0,t.max]
	h.hat.t_full <- as.matrix(h.hat.t$out_haz_full)
	# h.hat at observed unique time t = T.grid
	h.hat.t <- as.matrix(h.hat.t$out_haz)
	# ================================================================================================
	# estimate censoring G(A|W)
	# ================================================================================================
	message('estimating censoring')
	G.hat.t <- censor_SL_wrapper(dat = dat, T.uniq = T.uniq,
	                             Delta.SL.Lib = Delta.SL.Lib)

	Gn.A1.t_full <- as.matrix(G.hat.t$out_censor_full)
	Gn.A1.t <- as.matrix(G.hat.t$out_censor)
	# ================================================================================================
	# Gn.A1.t
	# ================================================================================================
	# plot initial fit
	if (verbose) lines(colMeans(Gn.A1.t) ~ T.uniq, col = 'yellow', lty = 1)
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
	# ================================================================================================
	# update
	# ================================================================================================
	stopping.criteria <- sqrt(l2.inner.step(Pn.D1.t, Pn.D1.t, T.uniq))/length(T.uniq) # 10-17

	update.tensor <- matrix(0, nrow = n.data, ncol = length(T.uniq))
	iter.count <- 0
	stopping.prev <- Inf

	# while ((stopping.criteria >= tol) & (iter.count <= max.iter)) { # ORGINAL
    while ((stopping.criteria >= tol) & (iter.count <= max.iter) & ((stopping.prev - stopping.criteria) >= max(-tol, -1e-5))) { #WILSON: TEMPORARY
		print(stopping.criteria)
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
		qn.current_full <- qn.A1.t_full * exp(epsilon.step * replicate(T.max, intergrand[,1])) #10-23

		# normalize the updated qn
		norm.factor <- compute.step.cdf(pdf.mat = qn.current, t.vec = T.uniq, start = Inf)[,1] #09-06
		qn.current[norm.factor > 1,] <- qn.current[norm.factor > 1,] / norm.factor[norm.factor > 1] #09-06
		qn.current_full[norm.factor > 1,] <- qn.current_full[norm.factor > 1,] / norm.factor[norm.factor > 1] #10-23

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
		# =============================================================================
		# compute new h_t
		h.hat.t_full_current <- matrix(0, nrow = n.data, ncol = max(T.uniq))
		for (it.n in 1:n.data) {
		    h.hat.t_full_current[it.n, ] <- qn.current_full[it.n, ] / Qn.current_full[it.n,]
		}
		# ------------------------------------------------
		# print('New D1')
		# compute new D1
		updated_IC <- compute_IC(dat = dat,
		                         dW = dW,
		                         T.uniq = T.uniq,
		                         h.hat.t_full = h.hat.t_full_current,
		                         g.fitted = g.fitted,
		                         Gn.A1.t_full = Gn.A1.t_full,
		                         Qn.A1.t = Qn.current,
		                         Qn.A1.t_full = Qn.current_full)

		D1.t <- updated_IC$D1.t
		D1.A1.t <- updated_IC$D1.A1.t

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
	variables <- list(T.uniq)
	params <- list(stopping.criteria, epsilon.step, iter.count, max.iter)
	initial_fit <- list(h.hat.t, Qn.A1.t, qn.A1.t)
	to.return <- list(Psi.hat = Psi.hat,
	                  T.uniq = T.uniq,
	                  params = params,
	                  variables = variables,
	                  initial_fit = initial_fit)
	class(to.return) <- 'onestep.surv'
	return(to.return)

}
