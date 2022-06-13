# drda

[![Build Status](https://github.com/albertopessia/drda/actions/workflows/r-cmd-check.yml/badge.svg?branch=master)](https://github.com/albertopessia/drda/actions/workflows/r-cmd-check.yml) [![Coverage](https://codecov.io/gh/albertopessia/drda/branch/master/graph/badge.svg)](https://codecov.io/gh/albertopessia/drda)

## Overview

*drda* is a [R](https://www.r-project.org/) package for fitting growth curves
and performing dose-response data analysis.

The current available models are:

- 5-parameter logistic function
- 5-parameter log-logistic function
- 4-parameter logistic function
- 4-parameter log-logistic function
- 2-parameter logistic function
- 2-parameter log-logistic function
- Gompertz function
- log-Gompertz function

A 6-parameter (log-)logistic function is also available for theoretical
research, but its use in real applications is discouraged because it is usually
non-identifiable from data.

## Installation

You can install the latest stable release of the `drda` R package from CRAN with

```{r}
install.packages("drda")
```

To install the latest development version of the package from GitHub use

```{r}
# install.packages("remotes")
remotes::install_github("albertopessia/drda")
```

## Usage

### Load the package

```{r}
library(drda)
```

### Example data

Package `drda` comes with our own example dataset:

```{r}
?voropm2

head(voropm2)
```

### Default fitting

```{r}
# by default `drda` uses a 4-parameter logistic function for model fitting

# common R API for fitting models
fit <- drda(response ~ log_dose, data = voropm2)

# get a general overview of the results
summary(fit)

# get parameter estimates by using generic functions...
coef(fit)
sigma(fit)

# ... or accessing the variables directly
fit$coefficients
fit$sigma

# compare the estimated model against a flat horizontal line using AIC, BIC, or
# the Likelihood Ratio Test (LRT)
#
# note that the LRT is testing the null hypothesis of a flat horizontal line
# being as a good fit as the full model, therefore we expect the test to be
# significant
#
# if the test is not significant, a horizontal line is probably a better model
anova(fit)
```

### Other models

```{r}
# use the `mean_function` argument to select a different model
fit_l2 <- drda(response ~ log_dose, data = voropm2, mean_function = "logistic2")
fit_l4 <- drda(response ~ log_dose, data = voropm2, mean_function = "logistic4")
fit_l5 <- drda(response ~ log_dose, data = voropm2, mean_function = "logistic5")
fit_gz <- drda(response ~ log_dose, data = voropm2, mean_function = "gompertz")

# which model should be chosen?
anova(fit_l2, fit_l4, fit_l5, fit_gz)

# 5-parameter logistic function provides the best fit (AIC and BIC are minimum)
```

### Weighted fit

```{r}
# it is possible to give each observation its own weight
fit_weighted <- drda(response ~ log_dose, data = voropm2, weights = weight)

# all the commands shown so far are available for a weighted fit as well
```

### Constrained optimization

```{r}
# it is possible to fix parameter values by setting the `lower_bound` and
# `upper_bound` appropriately
#
# unconstrained parameters have a lower bound of `-Inf` and an upper bound of
# `Inf`
#
# Important: be careful when deciding the constraints, because the optimization
#            problem might become very difficult to solve within a reasonable
#            number of iterations.
#
# In this particular example we are:
#   - fixing the alpha parameter to 1
#   - fixing the delta parameter to -1
#   - constraining the growth rate to be between 1 and 5
#   - not constraining the `phi` parameter, i.e. the `log(EC50)`
lb <- c(1, -1, 1, -Inf)
ub <- c(1, -1, 5,  Inf)

fit <- drda(
  response ~ log_dose, data = voropm2, lower_bound = lb, upper_bound = ub,
  max_iter = 260
)

summary(fit)

# if the algorithm does not converge, we can try to increase the maximum number
# of iterations or provide our own starting point
fit <- drda(
  response ~ log_dose, data = voropm2, lower_bound = lb, upper_bound = ub,
  start = c(1, -1, 2.6, 5), max_iter = 10000
)

summary(fit)

# our starting point was actually worse than the automatic one.
```

### Basic plot functionality

```{r}
fit_l5 <- drda(response ~ log_dose, data = voropm2, mean_function = "logistic5")

# plot the data used for fitting, the maximum likelihood curve, and
# *approximate* confidence intervals for the curve
plot(fit_l5)

# combine all curves in the same plot
fit_l2 <- drda(response ~ log_dose, data = voropm2, mean_function = "logistic2")
fit_l4 <- drda(response ~ log_dose, data = voropm2, mean_function = "logistic4")
plot(fit_l2, fit_l4, fit_l5)

# modify default plotting options
# use `legend_show = FALSE` to remove the legend altogether
plot(
  fit_l5, base = "10", col = "magenta", xlab = "x", ylab = "y", level = 0.9,
  midpoint = FALSE, main = "Example plot", legend_location = "topright",
  legend = "5-parameter logistic function"
)
```

## License

This package is free and open source software licensed under [MIT](LICENSE).
