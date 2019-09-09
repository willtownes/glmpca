# R package: glmpca

[![Build Status](https://travis-ci.com/willtownes/glmpca.svg?token=o1x5ZKVR5sA6MpqhDnQX&branch=master)](https://travis-ci.com/willtownes/glmpca)
[![codecov](https://codecov.io/gh/willtownes/glmpca/branch/master/graph/badge.svg?token=0bpQ61gRFj)](https://codecov.io/gh/willtownes/glmpca)

Generalized PCA for non-normally distributed data. If you find this useful please cite [Feature Selection and Dimension Reduction based on a Multinomial Model](https://www.biorxiv.org/content/10.1101/574574v1).

## Installation

```r
devtools::install_github("willtownes/glmpca")
```

## Usage

```r
library(glmpca)
res<-glmpca(some_data_matrix_obs_in_cols, 2)
factors<-res$factors
plot(factors[,1],factors[,2])
```

For more details see the vignettes.

## Issues and bug reports

Please use https://github.com/willtownes/glmpca/issues to submit issues, bug reports, and comments.
