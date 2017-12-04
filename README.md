# onestep.survival

<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) -->
<!-- [![](http://cranlogs.r-pkg.org/badges/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) -->
<!-- [![Travis-CI Build Status](https://travis-ci.org/wilsoncai1992/onestep.survival.svg?branch=master)](https://travis-ci.org/wilsoncai1992/onestep.survival) -->

<img style="float: left;margin:0 5rem 0 0" src="http://media.web.britannica.com/eb-media/29/76829-050-CD9C4B43.jpg" width="30%" height="30%">
<br>
<!-- <img style="float: left;margin:0 5rem 0 0" src="http://www.feenixx.com/space-exploration/posters/First_Step_on_Moon_Poster.jpg" width="30%" height="30%">
<br>
 -->
The `onestep.survival` R package is a tool for estimating counterfactual survival curve under static or dynamic interventions on treatment (exposure), while at the same time adjust for *measured* counfounding. Targeted Maximum Likelihood Estimate (TMLE) approach is employed to create a doubly robust and semi-parametrically efficient estimator. Machine Learning algorithms (SuperLearner) are implemented to all stages of the estimation.

Currently implemented **estimators** include:

1. One-step TMLE for the whole survival curve
2. One-step TMLE for survival at a specific end point
3. Iterative TMLE for survival at a specific end point

## Installation

To install the development version (requires the `devtools` package):

```R
install.packages('SuperLearner')
install.packages('tmle')
install.packages('Matrix')
devtools::install_github('wilsoncai1992/survtmle2')
devtools::install_github('wilsoncai1992/onestep.survival')
```

## Documentation

* To see all available package documentation:

```R
?onestep.survival
help(package = 'onestep.survival')
```

## Brief overview

### Data structure

The data input of all methods in the package should be an `R` `data.frame` in the following survival long data format:

```R
#   ID W A T.tilde delta
# 1  1 0 0      95     1
# 2  2 1 1       1     0
# 3  3 0 0     215     1
# 4  4 1 1      15     1
# 5  5 0 0      73     1
# 6  6 0 0      15     1
```

### Usage

```R
# simulate data
library(simcausal)
D <- DAG.empty()
D <- D +
  node("W", distr = "rbinom", size = 1, prob = .5) +
  node("A", distr = "rbinom", size = 1, prob = .35 + .4*W) +
  # Time to failure has an exponential distribution:
  node("Trexp", distr = "rexp", rate = 1 + .5*W - .5*A) +
  # Time to censoring has a weibull distribution:
  node("Cweib", distr = "rweibull", shape = .7 - .2*W, scale = 1) +
  # Actual time to failure is scaled by 100:
  node("T", distr = "rconst", const = round(Trexp*100,0)) +
  node("C", distr = "rconst", const = round(Cweib*100, 0)) +
  # Observed random variable (follow-up time):
  node("T.tilde", distr = "rconst", const = ifelse(T <= C , T, C)) +
  # Observed random variable (censoring indicator, 1 - failure event, 0 - censored):
  node("delta", distr = "rconst", const = ifelse(T <= C , 1, 0))
setD <- set.DAG(D)

# Simulate the data from the above data generating distribution:
dat <- sim(setD, n=1e2)
# subset into observed dataset
library(dplyr)
# only grab ID, W's, A, T.tilde, Delta
Wname <- grep('W', colnames(dat), value = TRUE)
dat <- dat[,c('ID', Wname, 'A', "T.tilde", "delta")]
head(dat)
# check positivity
check_positivity(dat)

library(onestep.survival)
# iterative TMLE: each time separately
dW <- rep(1, nrow(dat))
# dW <- rep(0, nrow(dat))
iterative_tmle <- survtmle_multi_t(dat = dat, dW = dW,
                                  SL.ftime = c("SL.glm","SL.mean","SL.step"),
                                  SL.ctime = c("SL.glm","SL.mean"),
                                  SL.trt = c("SL.glm","SL.mean","SL.step"))
plot(iterative_tmle, add = TRUE)

# one-step TMLE: target entire curve
dW <- rep(1, nrow(dat))
# dW <- rep(0, nrow(dat))
onestep_whole_curve <- surv_onestep(dat = dat, dW = dW,
                                    g.SL.Lib = c("SL.glm","SL.mean","SL.step"),
                                    Delta.SL.Lib = c("SL.glm","SL.mean", 'SL.gam'),
                                    ht.SL.Lib = c("SL.glm","SL.mean","SL.step", 'SL.gam'))
plot(onestep_whole_curve, add = TRUE)
```



## Citation
To cite `onestep.survival` in publications, please use:
> Cai W, van der Laan MJ (2016). *One-step TMLE for time-to-event outcomes.* Working paper.

## Funding

## Copyright
This software is distributed under the GPL-2 license.

## Community Guidelines
Feedback, bug reports (and fixes!), and feature requests are welcome; file issues or seek support [here](https://github.com/wilsoncai1992/onestep_survival/issues).