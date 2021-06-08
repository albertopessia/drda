# cran-comments

## drda 1.0.0

This is a resubmission for a new package.

- Added citation to our paper in `inst/CITATION` and `DESCRIPTION`.
- Added the manuscript as a vignette.
- Edited the `LICENSE` file to follow the CRAN template.

## Test environments
- local Windows 10 on R 4.1.0
- Ubuntu 16.04.6 LTS on Travis CI (r-oldrel, r-release, r-devel)
- win-builder x86_64-w64-mingw32 (R unstable 2021-06-07 r80458)
- R-hub fedora-clang-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub windows-x86_64-release (r-release)

## R CMD check results
> * checking CRAN incoming feasibility ... NOTE
> Maintainer: 'Alberto Pessia <dev@albertopessia.com>'
>
> New submission
>
> Possibly mis-spelled words in DESCRIPTION:
>   Malyutina (21:33)
>   Pessia (21:60)

The possibly mis-spelled words are author names.

> On windows-x86_64-release (r-release)
> checking sizes of PDF files under 'inst/doc' ... NOTE
> Unable to find GhostScript executable to run checks on size reduction

This note refers to the R-hub Windows platform and is not package-related.
The note does not appear in all other testing environments.
