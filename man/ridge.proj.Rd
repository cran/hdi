\name{ridge.proj}
\alias{ridge.proj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{P-values based on ridge projection method}
\description{P-values based on ridge projection method}
\usage{
ridge.proj(x, y, family = "gaussian", standardize = TRUE,
           lambda = 1, sigma = NULL, multiplecorr.method = "holm",
           N = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{family}{family}
  \item{standardize}{Should design matrix be standardized to unit column
    standard deviation (logical)?}
  \item{lambda}{Value of penalty parameter lambda (ridge regression).}
  \item{sigma}{Estimate of error standard deviation}
  \item{multiplecorr.method}{Either "WY" or any of
    \code{\link{p.adjust.methods}}.}
  \item{N}{Number of empirical samples (only used if multiplecorr.method
    == "WY")}
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

\references{B\"uhlmann, P. (2013), \emph{Statistical significance in
    high-dimensional linear models}, Bernoulli 19, 1212-1242.}

\author{Peter Buehlmann, Ruben Dezeure, Lukas Meier}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
x <- matrix(rnorm(100*100), nrow = 100, ncol = 100)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
fit.ridge <- ridge.proj(x, y)
which(fit.ridge$pval.corr < 0.05)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
