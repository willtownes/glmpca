# R package: glmpca

[![Build Status](https://travis-ci.com/willtownes/glmpca.svg?token=o1x5ZKVR5sA6MpqhDnQX&branch=master)](https://travis-ci.com/willtownes/glmpca)
[![codecov](https://codecov.io/gh/willtownes/glmpca/branch/master/graph/badge.svg)](https://codecov.io/gh/willtownes/glmpca)

Generalized PCA for non-normally distributed data. If you find this useful please cite [Feature Selection and Dimension Reduction based on a Multinomial Model.](https://doi.org/10.1186/s13059-019-1861-6) (doi:10.1186/s13059-019-1861-6)

A [python implementation](https://github.com/willtownes/glmpca-py) is also available.

## Installation

The [glmpca package](https://CRAN.R-project.org/package=glmpca) is available from CRAN. To install the stable release (recommended):

```r
install.packages("glmpca")
```

To install the development version:

```r
remotes::install_github("willtownes/glmpca")
```

## Usage

```r
library(glmpca)

#create a simple dataset with two clusters
mu<-rep(c(.5,3),each=10)
mu<-matrix(exp(rnorm(100*20)),nrow=100)
mu[,1:10]<-mu[,1:10]*exp(rnorm(100))
clust<-rep(c("red","black"),each=10)
Y<-matrix(rpois(prod(dim(mu)),mu),nrow=nrow(mu))

#visualize the latent structure
res<-glmpca(Y, 2)
factors<-res$factors
plot(factors[,1],factors[,2],col=clust,pch=19)
```

For more details see the vignettes. For compatibility with Bioconductor, see
[scry](https://bioconductor.org/packages/release/bioc/html/scry.html). 
For compatibility with Seurat objects, 
see [Seurat-wrappers](https://github.com/satijalab/seurat-wrappers).

## Alternative implementations

GLM-PCA has been around for awhile and we have not been able to dedicate as much time to its maintenance and ongoing improvement as we would like. Fortunately, there are numerous alternative implementations that improve on our basic idea. Many of them are likely to be faster and more memory-efficient than our version, and some have interesting additional capabilities such as uncertainty quantification. In reverse chronological order, here are some packages to check out.

* [fastglmpca](https://github.com/stephenslab/fastglmpca). Preprint:  [Weine, Carbonetto, & Stephens (2024)](https://doi.org/10.1101/2024.03.23.586420).
* [scGBM](https://github.com/phillipnicol/scGBM). Preprint: [Nicol & Miller 2023](https://doi.org/10.1101/2023.04.21.537881).
* [NewWave](https://bioconductor.org/packages/release/bioc/html/NewWave.html). Publication: [Agostinis et al 2022](https://doi.org/10.1093/bioinformatics/btac149).
* [LDVAE](https://docs.scvi-tools.org/en/stable/user_guide/models/linearscvi.html). Publication: [Svensson et al 2020](https://doi.org/10.1093/bioinformatics/btaa169).


## Issues and bug reports

Please use https://github.com/willtownes/glmpca/issues to submit issues, bug reports, and comments.
