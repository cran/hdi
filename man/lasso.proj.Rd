\name{lasso.proj}
\alias{lasso.proj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{P-values based on lasso projection method}
\description{P-values based on lasso projection method}
\usage{
lasso.proj(x, y, family = "gaussian", standardize = TRUE,
           multiplecorr.method = "holm", N = 10000,
           parallel = FALSE, ncores = 4,
           sigma = NULL, Z = NULL)
}
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{family}{family}
  \item{standardize}{Should design matrix be standardized to unit column
    standard deviation.}
  \item{multiplecorr.method}{Either "WY" or any of
    \code{\link{p.adjust.methods}}.}
  \item{N}{Number of empirical samples (only used if multiplecorr.method
  == "WY")}
  \item{parallel}{Should parallelization be used? (logical)}
  \item{ncores}{Number of cores used for parallelization.}
  \item{sigma}{Estimate of standard deviation of error term.}
  \item{Z}{user input}
}
\value{
  \item{pval}{Individual p-values for each parameter.}
  \item{pval.corr}{Multiple testing corrected p-values for each
    parameter.}
  \item{groupTest}{Function to perform groupwise tests. Groups are
    indicated using an index vector with entries in {1,...,p}.}
  \item{clusterGroupTest}{Function to perform groupwise tests based on
    hierarchical clustering. You can either provide a distance matrix
    and clustering method or the output of hierarchical clustering from
    the function \code{\link{hclust}} as for \code{\link{clusterGroupBound}}.}
  %\item{betahat}{initial estimate by the scaled lasso of \eqn{\beta^0}}
  %\item{bhat}{de-sparsified \eqn{\beta^0} estimate used for p-value calculation}
  \item{sigmahat}{\eqn{\widehat{\sigma}} coming from the scaled lasso.}
}
\references{van de Geer, S., B\"uhlmann, P., Ritov, Y. and Dezeure,
  R. (2013). On asymptotically optimal confidence regions and tests for
  high-dimensional models. Preprint.

  Zhang, C.,  Zhang, S. (2014), Confidence intervals for low dimensional
  parameters in high dimensional linear models, Journal of the Royal
  Statistical Society: Series B (Statistical Methodology).}
\author{Ruben Dezeure}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
x <- matrix(rnorm(100*20), nrow = 100, ncol = 20)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
fit.lasso <- lasso.proj(x, y)
which(fit.lasso$pval.corr < 0.05)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
