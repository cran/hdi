\name{groupBound}
\alias{groupBound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Lower bound on the l1-norm of groups of regression variables}

\description{In a (high-dimensional) regression, the function returns a
 lower bound that forms a one-sided confidence interval for the group
 l1-norm of a specified group of regression parameters.  It is assumed
 that errors have a Gaussian distribution with unknown noise level. The
 underlying vector that inference is made about is the l1-sparsest
 approximation to the noiseless data. Under a weak compatbility
 condition, this is identical to inference about the l1-sparsest
 approximation to the noiseless data.}

\usage{
groupBound(x, y, group, alpha = 0.05, nsplit = 11,
                s = min(10, ncol(x) - 1), setseed = TRUE,
                silent = FALSE, lpSolve = TRUE,
                parallel = FALSE, ncores = 4)}
\arguments{
  \item{x}{The design matrix of the regression with p columns for p predictor
    variables and n rows that correspond to n observations.}
  \item{y}{The response variable; a numeric vector of length n.}
  \item{group}{Either a numeric vector with entries in {1,...,p} or a
    list with such numeric vectors. If \code{group} is just a numeric
    vector, this is the group of variables for which a lower bound is
    computed. If \code{group} is a list, the lower bound is computed for
    each group in the list.}
  \item{alpha}{The level at which the test/ confidence interval is
    computed; a numeric value in (0,1).}
  \item{nsplit}{The number of data splits used.}
  \item{s}{The dimensionaility of the projection that is used. Lower
    values lead to faster computation and if n>50, then s is set to 50
    if left unspecified to avoid lengthy computations.}
  \item{setseed}{If \code{setseed} is true (recommended), then the same
    random seeds are used for all groups, which makes the confidence
    intervals simulatenously valid over all groups of variables tested.}
  \item{silent}{Output is supressed if this option is set to true.}
  \item{lpSolve}{Only set to false if lpSolve is not working on the
    current machine (setting it to false will results in much slower
    computations; only use on small problems).}
  \item{parallel}{Should parallelization be used? (logical)}
  \item{ncores}{Number of cores used for parallelization.}}

\details{The data are split since the noise level is unknown. On the
  first part of the random split, a cross-validated lasso solution is
  computed, using the glmnet implementation. This estimator is used as
  an initial estimator on the second half of the data. Results at level
  alpha are aggregated over \code{nsplit} splits via the median of
  results at levels alpha/2.}

\value{If \code{group} is a single numeric vector, a scalar containg the lower
 bound for this group of variables is returned. If \code{group} is a
 list, a numeric vector is retuned where each entry corresponds to the
 group of variables defined in the same order in \code{group}.}

\references{Nicolai Meinshausen (2013) 
 Assumption-free confidence intervals for groups of variables in sparse
 high-dimensional regression. http://arxiv.org/abs/1309.3489}
\author{Nicolai Meinshausen meinshausen@stat.math.ethz.ch}
%\note{
%}

\seealso{Use \code{clusterGroupBound} to test all groups in a hierarchical
 clustering tree.}

\examples{
## Create a regression problem with block-design: p = 10, n = 30,
## block size B = 5 and within-block correlation of rho = 0.99
p   <- 10
n   <- 100
B   <- 5
rho <- 0.99

ind <- rep(1:ceiling(p / B), each = B)[1:p]
Sigma <- diag(p)

for (ii in unique(ind)){
  id <- which(ind == ii)
  Sigma[id, id] <- rho
}
diag(Sigma) <- 1

x <- matrix(rnorm(n * p), nrow = n) \%*\% chol(Sigma)

## Create response with active variables 1 and 21 
beta    <- rep(0, p)
beta[1] <- 5

y  <- as.numeric(x \%*\% beta + rnorm(n))
            
## Compute lower bounds:

## Lower bound for the l1-norm of all variables 1-10 of the sparsest
## optimal vector  
lowerBoundAll <- groupBound(x, y, 1:p)
print(lowerBoundAll)
cat("\nlower bound for all variables 1-10: ", lowerBoundAll, "\n")

## Compute lower bounds:
## Lower bounds for variable 1 in itself, then groups 1-5
lowerBound <- groupBound(x, y, list(1, 1:5))
cat("lower bound for the groups {1}, {1,...,5}: ", lowerBound, "\n")
}
\keyword{confidence intervals}
\keyword{regression}
