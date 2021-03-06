% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rmvn_bcm}
\alias{rmvn_bcm}
\alias{rmvn_rue}
\title{Efficient Sampling of Normal Posterior Distributions}
\usage{
rmvn_bcm(Phi, d, alpha)

rmvn_rue(Phi, d, alpha)
}
\arguments{
\item{Phi}{the \code{n x p} matrix \eqn{\Phi} from the posterior}

\item{d}{the \code{p}-dimensional vector \eqn{d=diag(D)}}

\item{alpha}{the \code{n}-dimensional vector from the posterior}
}
\value{
A \code{p x 1} matrix with a draw from the posterior distribution
}
\description{
The \code{rmvn_bcm} and \code{rmvn_rue} functions allow for efficient 
sampling from normal posterior distributions. The posterior distribution
should have the form \eqn{\mu = \Sigma\Phi'\alpha}, where the posterior
covariance matrix is \eqn{\Sigma = (\Phi'\Phi+D^{-1})^{-1}} and \eqn{D} is 
a diagonal matrix.
}
\details{
The \code{rmvn_bcm} is appropriate when \code{n<p}, whereas \code{rmvn_rue} 
is typically the faster alternative when \code{n>p}.

The sampling routines are based on the proposals by Bhattacharya, Chakraborty 
and Mallick (2016) and Rue (2001). The former is based on 
an idea of avoiding operations in the \eqn{p} dimension in favor of working
in the \eqn{n} dimension, which is why it is preferrable when \eqn{n<p}. The
latter sampling routine instead does the opposite.

The C++ implementations are available as headers and can therefore be called
directly in C++ (e.g. via Rcpp) if necessary by other packages.

The functions can be used for standard Bayesian linear regression. Let the
likelihood be \deqn{y|X, \beta, \sigma^2\sim N(X\beta, \sigma^2 I_n),} suppose
\eqn{\sigma^2} is known and let the prior be \deqn{\beta_j|\sigma^2\sim N(0, \lambda_j^2\sigma^2),}
where \eqn{\lambda_j^2} is a known constant. Let \eqn{\Lambda} be a diagonal matrix
with \eqn{\lambda_j^2} along the diagonal. To sample from the posterior
of \eqn{\beta}, we let \eqn{\Phi = X/\sigma}, \eqn{d = \sigma^2\diag(Lambda)}
and \eqn{\alpha = y/\sigma}. For more elaborate priors (such as scale mixtures),
the \eqn{\lambda_j^2} are no longer fixed constants but sampled as well 
in Markov Chain Monte Carlo. The two functions in the package can in that case
be used as Gibbs sampling steps to sample from the conditional posterior distribution
\eqn{\beta|X, y, \sigma^2, \Lambda}.
}
\examples{
X <- matrix(rnorm(1000), 50, 20)
d <- runif(20, 0, 1)
alpha <- rnorm(50)

rue <- rmvn_rue(X, d, alpha)
bcm <- rmvn_bcm(X, d, alpha)
}
\references{
Bhattacharya, A., Chakraborty, A. and Mallick, B. (2016) Fast 
sampling with Gaussian scale mixture priors in high-dimensional regression,
\emph{Biometrika}, 103(4):985-991, doi:10.1093/biomet/asw042

Rue, H. (2001) Fast sampling of Gaussian Markov random fiels, \emph{Journal
of the Royal Statistical Society: Series B}, 63, 325-339, doi:10.1111/1467-9868.00288
}
