# cran-comments

## 2022-07-06

- This is a re-submission for fixing errors at <https://cran.r-project.org/web/checks/check_results_drda.html>.
- Increased version to 2.0.1 because submission would not be accepted otherwise.

## 2022-06-16

- Initiating submission process

## drda 2.0.0

This release extends the available models by implementing the log-logistic
family (on top of the previously available logistic family).

New features also included in this release:

- New model parameterization, hence the decision to increase the major version
- Improved initialization algorithm
- Exported functions for evaluating gradient and Hessian
- New function `effective_dose`
- Added examples in the docs

## Test environments

- Windows 10, local, R 4.2.0
- macOS-11, GitHub Actions, R 4.2.0
- Windows Server 2022, GitHub Actions, R 4.2.0
- Ubuntu 20.04, GitHub Actions, R 4.1.3, R 4.2.0, r-devel

## R CMD check results

There were no NOTEs, ERRORs or WARNINGs.

## Reverse dependencies

`devtools::revdep()` showed that `drda` is not currently a dependency of any
other package.
