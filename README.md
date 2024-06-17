# HGWR: Hierarchical and Geographically Weighted Regression

This is an C++ implementation of Hierarchical and Geographically Weighted Regression (HGWR) model.
HGWR model divides coefficients into three types: local fixed effects, global fixed effects, and random effects.
If data have spatial hierarchical structures (especially are overlapping on some locations), it is worth trying this model to reach better fitness.
For more information about this model, please turn to related articles.

## Dependency

Toolchain:

- CMake (>= 3.12.0)

C++ dependency:

- [armadillo](http://arma.sourceforge.net/)
- [gsl](https://www.gnu.org/software/gsl/)


## Building

The whole package is managed by CMake.

Run the following scripts to configure and build:

```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release --target hgwrbml
```

Currently there is no installation script.

## Usage

Please refer Wiki page for further information.

## Related Articles

Hu, Yigong, Lu, Binbin, Ge, Yong, Dong, Guanpeng, 2022.
Uncovering spatial heterogeneity in real estate prices via combined hierarchical linear model and geographically weighted regression.
Environment and Planning B: Urban Analytics and City Science.
[DOI](https://journals.sagepub.com/doi/10.1177/23998083211063885)
