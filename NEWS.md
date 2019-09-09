# glmpca 0.1

Initial release of the GLM-PCA method of dimension reduction. The primary focus is on count data (Poisson and negative binomial likelihoods). Bernoulli likelihood has also been implemented for binary data but this has not been tested as thoroughly. Other features included in this release

* basic checking for error conditions such as the data matrix having rows or columns that are all zero
* L2 (ridge) penalty on latent variables to improve numerical stability
* Arbitrary row and/or column covariates may be included to remove nuisance variation from latent factors.
* negative binomial shape parameter automatically estimated.
* latent factors automatically projected, scaled, and rotated to remove correlation with covariates and maintain orthogonality and improve interpretability.