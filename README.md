# HGWR: Hierarchical and Geographically Weighted Regression

[![R Package Check](https://github.com/HPDell/hlmgwr-backfitting-ml/actions/workflows/R.yml/badge.svg?branch=master)](https://github.com/HPDell/hlmgwr-backfitting-ml/actions/workflows/R.yml)
[![CRAN](https://www.r-pkg.org/badges/version/hgwrr)](https://cran.r-project.org/package=hgwrr)

This is an C++ implementation of Hierarchical and Geographically Weighted Regression (HGWR) model.
HGWR model divides coefficients into three types: local fixed effects, global fixed effects, and random effects.
If data have spatial hierarchical structures (especially are overlapping on some locations), it is worth trying this model to reach better fitness.
For more information about this model, please turn to related articles.
An related R package is also provided.

## Dependency

Toolchain:

- CMake (>= 3.12.0)

C++ dependency:

- [armadillo](http://arma.sourceforge.net/) (>= 9.900)
- [gsl](https://www.gnu.org/software/gsl/)

R dependency:

- [Rcpp](https://cran.r-project.org/web/packages/Rcpp/) (>= 1.0.8)
- [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/)

Note that on Windows, libraries [local323.zip](https://www.stats.ox.ac.uk/pub/Rtools/goodies/multilib/local323.zip) need to be installed to R home folder.
And R tools is required to build and install the R package from source.

## Building

The whole package (including R package) is managed by CMake.

### C++ library

Run the following scripts to configure and build:

```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release --target hgwrbml
```

Currently there is no installation script.

### R package

Set `WITH_R=ON` in CMake cache and re-configure CMake project or run the following scripts:

```bash
cmake .. -DWITH_R=ON
```

Then in your CMake build folder, there will be a folder `hgwrr` which contains an R source package.
If you are going to install this package from source, just run:

```bash
R CMD INSTALL hgwrr
```

Then just load this package in your R script like:

```r
library(hgwrr)
```

Currently, this package is not available on CRAN. 

### Usage

Please refer Wiki page for further information.

## Related Articles

Hu, Yigong, Lu, Binbin, Ge, Yong, Dong, Guanpeng, 2022.
Uncovering spatial heterogeneity in real estate prices via combined hierarchical linear model and geographically weighted regression.
Environment and Planning B: Urban Analytics and City Science.
[DOI](https://journals.sagepub.com/doi/10.1177/23998083211063885)
