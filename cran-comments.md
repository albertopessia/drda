# cran-comments

## 2026-09-16

- Initiating submission process

## drda 2.1.0

This release removes the automatic fit of a 5-parameter (log-)logistic function
from `anova` and fixes numerical errors and bugs.

## Test environments

- Windows 11 25H2 (26200.9168), local, R 4.6.1
- macOS 26.6.2 (25G83), GitHub Actions, R 4.6.1
- Windows Server 2025 (10.0.26100), GitHub Actions, R 4.6.1
- Ubuntu 24.04.5 LTS, GitHub Actions, R 4.5.3, R 4.6.1, r-devel

## R CMD check results

There were no NOTEs, ERRORs or WARNINGs.

## Reverse dependencies

`devtools::revdep()` showed that `drda` is not currently a dependency of any
other package.
