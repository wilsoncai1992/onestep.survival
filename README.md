# onestep.survival

<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) -->
<!-- [![](http://cranlogs.r-pkg.org/badges/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) [![](http://cranlogs.r-pkg.org/badges/grand-total/onestep.survival)](http://cran.rstudio.com/web/packages/onestep.survival/index.html) -->
<!-- [![Travis-CI Build Status](https://travis-ci.org/wilsoncai1992/onestep.survival.svg?branch=master)](https://travis-ci.org/wilsoncai1992/onestep.survival) -->

The `onestep.survival` R package is a tool for estimating counterfactual survival curve under static or dynamic interventions on treatment (exposure), while at the same time adjust for *measured* counfounding. Targeted Maximum Likelihood Estimate (TMLE) approach is employed to create a doubly robust and semi-parametrically efficient estimator. Machine Learning algorithms (SuperLearner) are implemented to all stages of the estimation.

Currently implemented **estimators** include:

1. One-step TMLE for the whole survival curve
2. One-step TMLE for survival at a specific end point
3. Iterative TMLE for survival at a specific end point

## Installation

To install the CRAN release version of `onestep.survival`: 

```R
devtools::install_github('wilsoncai1992/survtmle')
install.packages('onestep.survival')
```

To install the development version (requires the `devtools` package):

```R
devtools::install_github('wilsoncai1992/survtmle')
devtools::install_github('wilsoncai1992/onestep.survival')
```

## Documentation

Once the package is installed, see the [vignette](https://cran.r-project.org/web/packages/onestep.survival/vignettes/onestep.survival_vignette.pdf), consult the internal package documentation and examples. 

* To see the vignette in R:

```R
vignette("onestep.survival_vignette", package="onestep.survival")
```

* To see all available package documentation:

```R
?onestep.survival
help(package = 'onestep.survival')
```

## Brief overview

### Data structure

The data input of all methods in the package should be an `R` `data.frame` in the following survival long data format:

```R
#   ID W1 W A   T delta
# 1  1  1 0 0  23     0
# 2  2  1 0 0  80     1
# 3  3  0 0 0  24     1
# 4  4  0 0 0 115     0
# 5  5  0 1 1   3     0
# 6  6  0 0 0   6     1
```

## Citation
To cite `onestep.survival` in publications, please use:
> Cai W, van der Laan MJ (2016). *onestep.survival: One-step TMLE for time-to-event outcomes.* R package version 0.1.

## Funding

## Copyright
This software is distributed under the GPL-2 license.

## Community Guidelines
Feedback, bug reports (and fixes!), and feature requests are welcome; file issues or seek support [here](https://github.com/wilsoncai1992/onestep_survival/issues).