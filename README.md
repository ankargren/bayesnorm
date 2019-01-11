
bayesnorm
=========

Efficient Sampling of Normal Posterior Distributions
----------------------------------------------------

The `rmvn_bcm` and `rmvn_rue` functions allow for efficient sampling from normal posterior distributions. The posterior distribution should have the form `mu= Sigma * Phi' * alpha`, where the posterior covariance matrix is `Sigma = (Phi' * Phi + D^{-1})^{-1}` and `D` is a diagonal matrix.

The `rmvn_bcm` is appropriate when `n<p`, whereas `rmvn_rue` is typically the faster alternative when `n>p`.

The sampling routines are based on the proposals by Bhattacharya, Chakraborty and Mallick (2016) and Rue (2001). The former is based on an idea of avoiding operations in the `p` dimension in favor of working in the `n` dimension, which is why it is preferrable when `n<p`. The latter sampling routine instead does the opposite.

The C++ implementations are available as headers and can therefore be called directly in C++ (e.g. via Rcpp) if necessary by other packages.

Installation
------------

Example
-------

Incorporation in other packages
-------------------------------

### References

Bhattacharya, A., Chakraborty, A. and Mallick, B. (2016) Fast sampling with Gaussian scale mixture priors in high-dimensional regression, *Biometrika*, 103(4):985-991, [doi:10.1093/biomet/asw042](https://doi.org/10.1093/biomet/asw042)

Rue, H. (2001) Fast sampling of Gaussian Markov random fiels, *Journal of the Royal Statistical Society: Series B*, 63, 325-339, [doi:10.1111/1467-9868.00288](https://doi.org/10.1111/1467-9868.00288)
