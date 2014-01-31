\name{plot.clusterLowerBound}
\alias{plot.clusterLowerBound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot output of hierarchical testing of groups of variables}
\description{The functions plots the outcome of applying a lower bound on the
 l1-norm on groups of variables in a hierarchical clustering tree.} 
\usage{\method{plot}{clusterLowerBound}(x, cexfactor = 1, yaxis = "members",
     col = NULL, ...)}
\arguments{
  \item{x}{The output of function \code{clusterLowerBound}}
  \item{cexfactor}{Multiplies the size of the node symbols.}
  \item{yaxis}{For the default value ("members"), the hierarchical tree
    is shown as function of cluster size on the y-axis, whereas the node
    sizes are proportional to the lower l1-norm of the respective groups
    of variables. If \code{yaxis} takes any different value, then this
    is reversed and the tree is shown against the lower l1-norm on the
    y-axis, while node sizes are now proportional to the number of
    elements in each cluster.}
  \item{col}{The colour of the symbols for the nodes.}
  \item{...}{Additional arguments.}
}
%\details{}
\value{Nothing is returned}
\references{
 Nicolai Meinshausen (2013) 
 Assumption-free confidence intervals for groups of variables in sparse
 high-dimensional regression. 
 http://arxiv.org/abs/1309.3489
}
\author{Nicolai Meinshausen meinshausen@stat.math.ethz.ch}
%\note{}

\seealso{
  Use \code{clusterLowerBound} to test all groups in a hierarchical
  clustering tree. 
  Use \code{lowerGroupBound} to compute the lower bound for selected
  groups of variables. 
}
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
            
## Compute the lower bound for all groups in a hierarchical clustering tree
out <- clusterLowerBound(x, y, nsplit = 5)

## Plot the tree with y-axis proportional to the (log) of the number of
## group members and node sizes proportional to the lower l1-norm bound.
plot(out)

## Show the lower bound on the y-axis and node sizes proportional to
## number of group members
plot(out, yaxis = "")
}
\keyword{confidence intervals}
\keyword{regression}
\keyword{hierarchical clustering}
