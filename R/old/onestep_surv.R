source('./eval.l2.fast.R')
# source('./eval.l2.step.R')
# epsilon.step <- 1e-3
epsilon.step <- 1e-2
max.iter <- 1e2
n.data <- nrow(dat)
tol <- 1/n.data

dW <- rep(1, nrow(dat))

W <- dat$W
W <- as.data.frame(W)
# ================================================================================================
library(SuperLearner)
SL.library<- c("SL.glm", "SL.step", "SL.glm.interaction" )
gHatSL <- SuperLearner(Y=dat$A, X=W, SL.library=SL.library, family="binomial")
# g.hat for each observation
g.fitted <- gHatSL$SL.predict
# ================================================================================================
# conditional hazard (by SL)
# ================================================================================================
T.uniq <- sort(unique(dat$T))
h.hat.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))

for (it1 in 1:length(T.uniq)) {
	print(it1)
	T.it <- T.uniq[it1]
	I.t <- (dat$T == T.it) + 0
	# ----------------------------------------------------------------------------------------
	# SL for h.n
	# ----------------------------------------------------------------------------------------
	X.mat <- dat[,c('A', 'W')]
	# predict counterfactual
	X.mat.new <- X.mat
	X.mat.new[,1] <- 1
	####################################
	# Temporaray!
	# subset for conditional hazard
	# # is.at.risk <- (dat$T >= T.it)
	# # X.mat <- X.mat[is.at.risk,]
	# # I.t <- I.t[is.at.risk]
	####################################
	# hHat.SL <- SuperLearner(Y=I.t, X=X.mat, newX = X.mat.new,
	# 												SL.library=SL.library, family="binomial")
	# hn.A1 <- hHat.SL$SL.predict
	# ----------------------------------------------------------------------------------------
	# use glm for faster
	X.mat <- cbind(I.t, X.mat)
	hhat.glm <- glm(I.t ~ A + W, family = 'binomial', data = X.mat)
	hn.A1 <- predict(hhat.glm, newdata = X.mat.new, type = 'response')
	# ----------------------------------------------------------------------------------------
	h.hat.t[,it1] <- hn.A1
}
# ================================================================================================
# Qn.A1.t
# ================================================================================================
source('./gen.step.func.R')
source('./compute.step.cdf.R')
Qn.A1.t <- matrix(0, nrow = n.data, ncol = length(T.uniq))
for (it.n in 1:n.data) {
	# print(it.n)
	cum.haz.vec <- rep(NA, length(T.uniq))
	# ------------------------------------------------------------------------------------
	# compute cumulative hazard
	# ------------------------------------------------------------------------------------
	# cumsum approach
	# cum.haz.vec <- compute.step.cdf(pdf.vec = h.hat.t[it.n,], t.vec = T.uniq, start = -Inf)
	# ------------------------------------------------------------------------------------
	# integrate approach
	cum.haz.vec <- cumsum(h.hat.t[it.n,]) #WILSON: WRONG!!!!!! integration

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
Pn.D1.t.func <- gen.step.func(y.vec = c(0, Pn.D1.t), t.vec = T.uniq)

D1.t.func <- gen.step.func(y.vec = c(1, D1.t[it.n,]), t.vec = T.uniq)
# curve(Qn.A1.t.func, from = 0, to = 10000)
# ================================================================================================
# update
# ================================================================================================
compute.update <- function(D1.t.func.prev, Pn.D1.func.prev) {
	# curve(Pn.D1.func.prev, from = 0, to = T.max)
	# curve(D1.t.func.prev, from = 0, to = T.max)
	
	# formula on p.30
	# result <- l2.inner(Pn.D1.func.prev, D1.t.func.prev, min = 0, max = T.max) /
							# sqrt(l2.inner(D1.t.func.prev, D1.t.func.prev, min = 0, max = T.max))
	# formula on p.28
	result <- l2.inner(Pn.D1.func.prev, D1.t.func.prev, min = 0, max = T.max) /
		sqrt(l2.inner(Pn.D1.func.prev, Pn.D1.func.prev, min = 0, max = T.max))
	
	if (is.na(result)) {
		result <- NA
	}
	return(result)
}
# ================================================================================================
T.max <- max(T.uniq)

stopping.criteria <- sqrt(l2.inner(Pn.D1.t.func, Pn.D1.t.func, min = 0, max = T.max))
update.seq <- matrix(0, nrow = n.data, ncol = 1)
iter.count <- 0
while ((stopping.criteria >= tol) & (iter.count <= max.iter)) {
	print(stopping.criteria)

	# update the qn
	update.vec <- rep(0, n.data)
	for (it.n in 1:n.data) {
		# print(it.n)
		D1.t.func <- gen.step.func(y.vec = c(1, D1.t[it.n,]), t.vec = T.uniq)
		update.val <- compute.update(D1.t.func.prev = D1.t.func,
																 Pn.D1.func.prev = Pn.D1.t.func)
		update.vec[it.n] <- update.val
	}
	update.seq <- cbind(update.seq, update.vec)

	intergrand <- rowSums(update.seq)
	intergrand[is.na(intergrand)] <- 0
	qn.current <- qn.A1.t * exp(epsilon.step * intergrand)

	# ------------------------------------------------
	# compute new Qn
	# Qn.current <- apply(qn.current, 1, function(x) compute.step.cdf(pdf.vec = x, t.vec = T.uniq, start = Inf))
	# Qn.current <- t(Qn.current)
	
	# check error
	# all.equal(compute.step.cdf(pdf.vec = qn.current[1,], t.vec = T.uniq, start = Inf), Qn.current[1,])
	# ------------------------------------------------
	Qn.current <- matrix(NA, nrow = n.data, ncol = length(T.uniq))
	for (it.n in 1:n.data) {
		Qn.current[it.n,] <- rev(cumsum( rev(qn.current[it.n,]) ))
	}
	# ------------------------------------------------
	
	
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
	Pn.D1.t.func <- gen.step.func(y.vec = c(0, Pn.D1.t), t.vec = T.uniq)
	
	D1.t.func <- gen.step.func(y.vec = c(1, D1.t[it.n,]), t.vec = T.uniq)
	# ================================================================================================
	# new stopping criteria
	stopping.criteria <- sqrt(l2.inner(Pn.D1.t.func, Pn.D1.t.func, min = 0, max = T.max))
	# interation count += 1
	iter.count <- iter.count + 1
}

# ================================================================================================
# compute the target parameter
# ================================================================================================
Psi.hat <- colMeans(Qn.current)

# ================================================================================================
#
# ================================================================================================
for (play in 1:n.data) {
	play <- play + 1
	plot(h.hat.t[play,] ~ T.uniq, type = 'l')
}

