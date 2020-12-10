# drda
[![Build Status](https://travis-ci.com/albertopessia/drda.svg?branch=master)](https://travis-ci.com/albertopessia/drda) [![Coverage](https://codecov.io/gh/albertopessia/drda/branch/master/graph/badge.svg)](https://codecov.io/gh/albertopessia/drda)

## Overview

*drda* is a [R](https://www.r-project.org/) package for analysing dose-response
data.

## Installation

The package is currently available only in its development version.

```{r}
# install.packages("devtools")
devtools::install_github("albertopessia/drda")
```

## Usage

```{r}
# reproduce this example
set.seed(1492019)

library(drda)

# set true parameter as `c(lower_bound, upper_bound, growth_rate, log_ec50)`
theta <- c(0.02, 0.86, -1, -2)
sigma <- 0.05

# our model expects doses on the log scale
log_doses <- log(rep(c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), each = 3))

observations <- rnorm(
  n = length(log_doses),
  mean = logistic4_function(log_doses, theta),
  sd = sigma
)

test_data <- data.frame(
  x = log_doses,
  y = observations
)

# `drda` uses the common R API for fitting models
fit <- drda(y ~ x, data = test_data)

# the previous call is also equivalent to `drda(observations ~ log_doses)`

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
# Note that LRT is testing the null hypothesis of a flat horizontal line being
# as a good fit as the full model, therefore we expect the test to be
# significant
#
# If the test is not significant, a horizontal line is probably a better model.
anova(fit)

# constrained optimization
#
# It is possible to fix parameter values by setting the `lower_bound` and
# `upper_bound` to the desired value.
# Unconstrained parameters have a lower bound of `-Inf` and an upper bound of
# `Inf`.
#
# Important: be careful when deciding the constraints, because the optimization
#            problem might become very difficult to solve within a reasonable
#            number of iterations.
#
# In this particular example we are:
#   - fixing the `lower_bound` to 0
#   - fixing the `upper_bound` to 1
#   - constraining the growth rate to be between -5 and 5
#   - not constraining the `log_ec50`
lb <- c(0, 1, -5, -Inf)
ub <- c(0, 1,  5,  Inf)

# this is a difficult problem
fit <- drda(
  y ~ x, data = test_data, lower_bound = lb, upper_bound = ub, max_iter = 100
)

# note that the algorithm did not converge within the required tolerance, but
# it is still very close to the solution
summary(fit)

# if the algorithm does not converge, you can try to increase the maximum number
# of iterations or provide your own starting point
fit <- drda(
  y ~ x, data = test_data, lower_bound = lb, upper_bound = ub,
  start = c(0, 1, -0.6, -2), max_iter = 50000
)

summary(fit)

# we now give weights to each observation and use a weighted least squares
# estimation
#
# the weights used here are simply `|y - m|^(-p)`, where `m` is the
# corresponding common mean and `p = 0.25`.
w <- unlist(by(observations, log_doses, function(z) abs(z - mean(z))^(-0.25)))

fit_logistic4 <- drda(y ~ x, data = test_data, weights = w)

summary(fit_logistic4)

# we now fit a 2-parameter logistic function by fixing the first two parameters
# and estimating the other two
fit_logistic2 <- drda(
  y ~ x, data = test_data, weights = w,
  lower_bound = c(0, 1, -Inf, -Inf),
  upper_bound = c(0, 1, Inf, Inf)
)

summary(fit_logistic2)

# which model should we choose?
anova(fit_logistic4, fit_logistic2)

# 4-parameter logistic function provides the best fit (AIC is minimum)
```

## License
This package is free and open source software licensed under [MIT](LICENSE).
