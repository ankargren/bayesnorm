
bayesnorm
=========

Efficient Sampling of Normal Posterior Distributions
----------------------------------------------------

The `rmvn_bcm` and `rmvn_rue` functions allow for efficient sampling from normal posterior distributions. The posterior distribution should have the form *μ* = *Σ**Φ*′*α*, where the posterior covariance matrix is *Σ* = (*Φ*′*Φ* + *D*<sup>−1</sup>)<sup>−1</sup> and *D* is a diagonal matrix.

The `rmvn_bcm` is appropriate when `n<p`, whereas `rmvn_rue` is typically the faster alternative when `n>p`.

The sampling routines are based on the proposals by Bhattacharya, Chakraborty and Mallick (2016) and Rue (2001). The former is based on an idea of avoiding operations in the *p* dimension in favor of working in the *n* dimension, which is why it is preferrable when *n* &lt; *p*. The latter sampling routine instead does the opposite.

The C++ implementations are available as headers and can therefore be called directly in C++ (e.g. via Rcpp) if necessary by other packages.

The functions can be used for standard Bayesian linear regression. Let the likelihood be
*y*|*X*, *β*, *σ*<sup>2</sup> ∼ *N*(*X**β*, *σ*<sup>2</sup>*I*<sub>*n*</sub>),
 suppose *σ*<sup>2</sup> is known and let the prior be
*β*<sub>*j*</sub>|*σ*<sup>2</sup> ∼ *N*(0, *λ*<sub>*j*</sub><sup>2</sup>*σ*<sup>2</sup>),
 where *λ*<sub>*j*</sub><sup>2</sup> is a known constant. Let $} be a diagonal matrix with *λ*<sub>*j*</sub><sup>2</sup> along the diagonal. To sample from the posterior of *β*, we let *Φ* = *X*/*σ*, $d = \\sigma^2\\diag(Lambda)$ and *α* = *y*/*σ*. For more elaborate priors (such as scale mixtures), the *λ*<sub>*j*</sub><sup>2</sup> are no longer fixed constants but sampled as well in Markov Chain Monte Carlo. The two functions in the package can in that case be used as Gibbs sampling steps to sample from the conditional posterior distribution *β*|*X*, *y*, *σ*<sup>2</sup>, *Λ*.

### References

Bhattacharya, A., Chakraborty, A. and Mallick, B. (2016) Fast sampling with Gaussian scale mixture priors in high-dimensional regression, *Biometrika*, 103(4):985-991, <doi:10.1093/biomet/asw042>

Rue, H. (2001) Fast sampling of Gaussian Markov random fiels, *Journal of the Royal Statistical Society: Series B*, 63, 325-339, <doi:10.1111/1467-9868.00288>
