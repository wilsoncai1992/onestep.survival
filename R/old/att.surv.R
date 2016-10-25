library(simcausal)
D <- DAG.empty()

D <- D +
	# node("W", distr = "rnorm", mean = 1, sd = .5) +
	# node("W", distr = "rpois", lambda=.1) +
	node("W", distr = "rbinom", size = 1, prob = .5) +
	node("A", distr = "rbinom", size = 1, prob = .5) +
	# Time to failure has an exponential distribution:
	node("Trexp", distr = "rexp", rate = 1 + W - .5*A) +
	# Actual time to failure is scaled by 100:
	node("T", distr = "rconst", const = round(Trexp*100,0))
setD <- set.DAG(D)

# Simulate the data from the above data generating distribution:
dat <- sim(setD, n=1e3)
head(dat)


dW <- 1

# ================================================================================
# truth
# ================================================================================
library(survival)
n.data <- nrow(dat)
km.fit <- survfit(Surv(T,rep(1, n.data)) ~ A, data = dat)
plot(km.fit)
haz.play <- epi.insthaz(km.fit, conf.level = 0.95)
plot(haz.play$time, haz.play$est, xlab = "Days",
		 ylab = "Instantaneous hazard", type = "b", pch = 16)


# ================================================================================