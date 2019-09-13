# R package: glmpca

[![Build Status](https://travis-ci.com/willtownes/glmpca.svg?token=o1x5ZKVR5sA6MpqhDnQX&branch=master)](https://travis-ci.com/willtownes/glmpca)
[![codecov](https://codecov.io/gh/willtownes/glmpca/branch/master/graph/badge.svg)](https://codecov.io/gh/willtownes/glmpca)

Generalized PCA for non-normally distributed data. If you find this useful please cite [Feature Selection and Dimension Reduction based on a Multinomial Model](https://www.biorxiv.org/content/10.1101/574574v1).

## Installation

```r
devtools::install_github("willtownes/glmpca")
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

For more details see the vignettes.

## Issues and bug reports

Please use https://github.com/willtownes/glmpca/issues to submit issues, bug reports, and comments.
