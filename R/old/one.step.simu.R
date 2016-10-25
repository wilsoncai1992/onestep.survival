# source('./onestep_surv_new_inner.R')
# load fitting code for one-step survival curve
source('./onestep_surv_new_inner_new_update.R')
# =============================================================================
# simulate
# =============================================================================
library(simcausal)
D <- DAG.empty()

D <- D +
	# node("W", distr = "rnorm", mean = 1, sd = .5) +
	node("W", distr = "rpois", lambda=.1) +
	# node("W", distr = "rbinom", size = 1, prob = .5) +
	node("A", distr = "rbinom", size = 1, prob = .5) +
	# Time to failure has an exponential distribution:
	# node("Trexp", distr = "rexp", rate = 1 + W - .5*A) +
	node("Trexp", distr = "rexp", rate = 1.5 - 1*A) +
	# Actual time to failure is scaled by 100:
	node("T", distr = "rconst", const = round(Trexp*100,0))
setD <- set.DAG(D)

# Simulate the data from the above data generating distribution:
# dat <- sim(setD, n=1e3)
dat <- sim(setD, n=50)
head(dat)

# =============================================================================
# truth
# =============================================================================
# q <- seq(0,10,.1)
# # truesurvExp <- 1 - pexp(q, rate = 1)
# truesurvExp <- 1 - pexp(q, rate = .5)
# plot(round(q*100,0), truesurvExp, type="l", cex=0.2, col = 'red')
# =============================================================================
# KM
# =============================================================================
# library(survival)
# n.data <- nrow(dat)
# km.fit <- survfit(Surv(T,rep(1, n.data)) ~ A, data = dat)


## library(rms)
## survplot(km.fit)

# lines(km.fit)
# =============================================================================
# estimation
# =============================================================================
# dW <- rep(0, nrow(dat))
dW <- rep(1, nrow(dat))

draw.fn <- function() {
	output <- surv.one.step(dat, dW)
}
# lines(output[[1]] ~ output[[2]], type = 'l', col = 'green')

library(animation)
saveGIF(draw.fn(), interval = 0.1)
# =============================================================================
# library(epiR)
# haz.play <- epi.insthaz(km.fit, conf.level = 0.95)
# plot(haz.play$time, haz.play$est, xlab = "Days",
# 		 ylab = "Instantaneous hazard", type = "b", pch = 16)
