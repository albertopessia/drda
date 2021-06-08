# drda

[![Build Status](https://travis-ci.com/albertopessia/drda.svg?branch=master)](https://travis-ci.com/albertopessia/drda) [![Coverage](https://codecov.io/gh/albertopessia/drda/branch/master/graph/badge.svg)](https://codecov.io/gh/albertopessia/drda)

## Overview

*drda* is a [R](https://www.r-project.org/) package for fitting growth curves
and performing dose-response data analysis.

The current available models are:

- 5-parameter logistic function
- 4-parameter logistic function
- 2-parameter logistic function
- Gompertz function

A 6-parameter logistic function is also available for theoretical research, but
its use in real applications is discouraged because it is usually
non-identifiable from data.

## Installation

Install the package from CRAN:

```{r}
install.packages("drda")
```

## Usage

### Example data

```{r}
dose <- rep(c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), each = 3)
log_dose <- log(dose)

relative_viability <- c(
  0.87736, 0.81284, 0.88311, 0.87349, 0.84577, 0.99942, 0.88896, 0.73554,
  0.84204, 0.51804, 0.51926, 0.50125, 0.25321, 0.08394, -0.00072, 0.04925,
  0.07080, 0.09143, 0.04110, -0.03615, 0.09256
)

weights <- c(
  1.10261, 1.01677, 1.06836, 0.96812, 0.81738, 1.11376, 1.08168, 0.99638,
  0.99288, 0.90101, 1.09023, 1.19479, 0.90741, 1.14378, 0.84185, 0.82132,
  0.85104, 0.99346, 0.98071, 1.10807, 0.94382
)

test_data <- data.frame(
  y = relative_viability,
  x = log_dose
)
```

### Load the package

```{r}
library(drda)
```

### Default fitting

```{r}
# by default `drda` uses a 4-parameter logistic function for model fitting

# common R API for fitting models (the following two statements are equivalent)
fit <- drda(relative_viability ~ log_dose)
fit <- drda(y ~ x, data = test_data)

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
fit_logistic2 <- drda(y ~ x, data = test_data, mean_function = "logistic2")
fit_logistic4 <- drda(y ~ x, data = test_data, mean_function = "logistic4")
fit_logistic5 <- drda(y ~ x, data = test_data, mean_function = "logistic5")
fit_gompertz <- drda(y ~ x, data = test_data, mean_function = "gompertz")

# which model should be chosen?
anova(fit_logistic2, fit_logistic4, fit_logistic5, fit_gompertz)

# 4-parameter logistic function provides the best fit (AIC is minimum)
#
# the true model was indeed a 4-parameter function with
# `theta = c(0.02, 0.86, -1, -2)` and `sigma = 0.05`
```

### Weighted fit

```{r}
# it is possible to give each observation its own weight
fit_weighted <- drda(y ~ x, data = test_data, weights = weights)

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
#   - fixing the `lower_bound` to 0
#   - fixing the `upper_bound` to 1
#   - constraining the growth rate to be between -5 and 5
#   - not constraining the `phi` parameter, i.e. the `log(EC50)`
lb <- c(0, 1, -5, -Inf)
ub <- c(0, 1,  5,  Inf)

fit <- drda(
  y ~ x, data = test_data, lower_bound = lb, upper_bound = ub, max_iter = 100
)

summary(fit)

# if the algorithm does not converge, we can try to increase the maximum number
# of iterations or provide our own starting point
fit <- drda(
  y ~ x, data = test_data, lower_bound = lb, upper_bound = ub,
  start = c(0, 1, -0.6, -2), max_iter = 50000
)

summary(fit)

# our starting point was actually worse than the automatic one.
```

### Basic plot functionality

```{r}
fit <- drda(y ~ x, data = test_data, mean_function = "logistic5")

# plot the data used for fitting, the maximum likelihood curve, and
# *approximate* confidence intervals for the curve
plot(fit)

# when the model is a 4-parameter logistic function, or a 2-parameter logistic
# function, the `phi` parameter is also shown
fit <- drda(y ~ x, data = test_data, mean_function = "logistic4")
plot(fit)

fit <- drda(y ~ x, data = test_data, mean_function = "logistic2")
plot(fit)
```

## License

This package is free and open source software licensed under [MIT](LICENSE).
