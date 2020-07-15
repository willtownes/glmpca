# News
All notable changes to glmpca will be documented in this file.
The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)

## [0.2.0] - 2020-07
This release should be backwards-compatible even though some arguments are now deprecated.
### Added
- AvaGrad optimizer which is faster and more numerically stable than Fisher scoring.
- Minibatch methods (memoized and stochastic gradients) to support large sparse Matrices
- Automatic tuning of optimization hyperparameters (learning rate or penalty)
- Negative binomial likelihood with feature-specific overdispersion parameters
- print and predict methods for fitted glmpca object
### Changed
- more detailed documentation of main `glmpca` function
- split off large blocks of code in the main `glmpca` function to separate functions
- organize results as an S3 class instead of a list
- moved arguments `penalty` and `verbose` into the `ctl` parameters list
- deprecated family `mult` and `bern`, replaced with family `binom`
- control parameter `eps` has now been renamed `tol`

## [0.1.0] - 2019-09
Initial release of the GLM-PCA method of dimension reduction. The primary focus is on count data (Poisson and negative binomial likelihoods). Bernoulli likelihood has also been implemented for binary data but this has not been tested as thoroughly. 
### Added
- basic checking for error conditions such as the data matrix having rows or columns that are all zero
- L2 (ridge) penalty on latent variables to improve numerical stability
- Arbitrary row and/or column covariates may be included to remove nuisance variation from latent factors.
- negative binomial shape parameter automatically estimated.
- latent factors automatically projected, scaled, and rotated to remove correlation with covariates and maintain orthogonality and improve interpretability.
