# drda 2.0.3

- Updated vignette after review from the Journal of Statistical Software.
- Fixed again the `x`-axis of the plot.

# drda 2.0.2

- Fixed a bug in the `plot` when doses were not previously sorted.
- Fixed a bug in the `x`-axis of the `plot` when data contained zeros.
- It is now possible to not plot data points by setting option `plot_data` to
`FALSE` (default `TRUE`).

# drda 2.0.1

- Small fixes to unit tests because they were not passing on specific systems.
- Updated vignette with new simulation results.

# drda 2.0.0

It is now possible to fit models using either the log-dose or the dose scale.

To accommodate this extension it was necessary to change the default model
parameterization, which now follows that of the Emax model
[(Macdougall, 2006)](https://doi.org/10.1007/0-387-33706-7_9).

Briefly, the 5-parameter logistic function is now defined as

```{r}
alpha + delta / (1 + nu * exp(-eta * (x - phi)))^(1 / nu)
```

Parameter `alpha` is the value of the function when `x` approaches `-Inf`.
Parameter `delta` is the (signed) height of the curve.
Parameter `eta > 0` represents the steepness (growth rate) of the curve.
Parameter `phi` is related to the mid-value of the function.
Parameter `nu` affects near which asymptote maximum growth occurs.

Similarly, the newly implemented log-logistic function (when `x >= 0`) is
defined as

```{r}
alpha + delta * (x^eta / (x^eta + nu * phi^eta))^(1 / nu)
```

Check the vignette (`vignette("drda", package = "drda")`) or the help page (`help(drda)`) to know more about the available models.

Here is a change log from previous version:

* Change parameterization to follow that of the Emax model.
* Implement the log-logistic family of models.
* Improve initialization algorithm to be more efficient and (hopefully) robust.
* Exported functions for evaluating theoretical gradient and Hessian of each
implemented model.
* Implement the `effective_dose` function for estimating effective doses.
* Added examples to help pages.
* Many minor bug fixes (too many to list them all).

# drda 1.0.0

First public release.
