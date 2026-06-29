# DiPMInd-AOS-2026-code

Implementations of distance-profile-based test statistics for assessing pairwise and mutual independence between random variables, developed in support of:

> Yaqing Chen and Paromita Dubey. "DiPMInd: Distance profile based mutual independence testing for random objects." Annals of Statistics (to appear). [(https://arxiv.org/abs/2412.06766)]

## Overview

This repository contains C++ (Rcpp) source code for computing nonparametric test statistics that compare an empirical joint distribution to the product of its marginals, using **n × n distance matrices** as input rather than raw observations. Two routines are provided:

| File | Function | Purpose |
|---|---|---|
| `teststatRcppsorted_fastest.cpp` | `teststatv2Rcpp_sorted2D(Dx, Dy)` | Bivariate (pairwise) independence test statistics |
| `mutualTeststatRcppFast.cpp` | `mutualTeststatRcpp_sorted3D(Dx, Dy, Dz)` | Trivariate (mutual/joint) independence test statistics |

Both functions return a length-8 numeric vector with four weighting variants, computed separately for the **diagonal** (matched-index, e.g. `j == k`) and **product** (all index combinations) terms:

```
jt_stat, jt_AD, jt_W, jt_ADsum, prod_stat, prod_AD, prod_W, prod_ADsum
```

- `stat` — unweighted squared deviation between the empirical joint CDF and the product of marginal CDFs
- `AD` — weighted by the product of marginal variance terms, `F(1-F)`
- `ADsum` — weighted by the sum of marginal variance terms
- `W` — weighted by the product of marginal CDFs

## Inputs

Both functions take **n × n distance (or kernel/closed-ball indicator) matrices** as `Rcpp::NumericMatrix`, one per variable:

- `Dx`, `Dy` — for the bivariate statistic
- `Dx`, `Dy`, `Dz` — for the trivariate statistic

Matrices are expected to be square and aligned by row/column index across variables (row `i` corresponds to the same observation in `Dx`, `Dy`, and `Dz`).

## Requirements

- R (≥ 3.6 recommended)
- [Rcpp](https://cran.r-project.org/package=Rcpp)
- A C++11-compatible compiler (e.g. via Rtools on Windows, Xcode command line tools on macOS, or build-essential on Linux)

## Usage

```r
library(Rcpp)

sourceCpp("src/teststatRcppsorted_fastest.cpp")
sourceCpp("src/mutualTeststatRcppFast.cpp")

# Dx, Dy, Dz: n x n distance matrices, one per variable
result2D <- teststatv2Rcpp_sorted2D(Dx, Dy)
result3D <- mutualTeststatRcpp_sorted3D(Dx, Dy, Dz)

names(result2D) <- c("jt_stat", "jt_AD", "jt_W", "jt_ADsum",
                      "prod_stat", "prod_AD", "prod_W", "prod_ADsum")
names(result3D) <- names(result2D)
```

See `examples/example_usage.R` for a runnable example with simulated data.

## Repository structure

```
DiPMInd-AOS-2026-code/
├── README.md
├── LICENSE
├── src/
│   ├── teststatRcppsorted_fastest.cpp
│   └── mutualTeststatRcppFast.cpp
├── R/
│   └── (optional R wrapper/helper functions)
├── examples/
│   └── example_usage.R
└── data/
    └── (optional example or simulated datasets)
```

## Citation

If you use this code, please cite:

```
Yaqing Chen and Paromita Dubey. "DiPMInd: Distance Profile based Mutual Independence testing for random objects." arXiv preprint arXiv:2412.06766 (Annals of Statistics (to appear)).
```

## License

Released under the MIT License — see `LICENSE` for details.

## Contact

Paromita Dubey — paromita@marshall.usc.edu
