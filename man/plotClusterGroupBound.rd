\name{plot.clusterGroupBound}
\title{Plot output of hierarchical testing of groups of variables}
\alias{plot.clusterGroupBound}
\concept{confidence intervals}
\concept{hierarchical clustering}
\description{
  The \code{\link{plot}()} method for \code{"\link{clusterGroupBound}"} objects
  plots the outcome of applying a lower bound on the l1-norm on groups of
  variables in a hierarchical clustering tree.
}
\usage{\method{plot}{clusterGroupBound}(x, cexfactor = 1, yaxis = "members",
     xlab = "", col = NULL, pch = 20, \dots)
}
\arguments{
  \item{x}{an object of \code{\link{class}} \code{"clusterGroupBound"},
    as resulting from \code{\link{clusterGroupBound}()}.}
  \item{cexfactor}{numeric expansion factor for the size of the node symbols.}
  \item{yaxis}{a string; for the default \code{"members"}, the hierarchical tree
    is shown as function of cluster size on the y-axis, whereas the node
    sizes are proportional to the lower l1-norm of the respective groups
    of variables. If \code{yaxis} takes any different value, then this
    is reversed and the tree is shown against the lower l1-norm on the
    y-axis, while node sizes are now proportional to the number of
    elements in each cluster.}
  \item{xlab}{label used for the x-axis; by default none.}
  \item{col}{the colour of the symbols for the nodes.}
  \item{pch}{the plot symbol (see \code{\link{points}}) of the symbols
    for the nodes.}
  \item{\dots}{optional additional arguments passed to
    \code{\link{plot.default}}.}
}
%\details{}
\value{Nothing is returned}
%\references{}
\author{Nicolai Meinshausen \email{meinshausen@stat.math.ethz.ch}}
%\note{}

\seealso{
  Use \code{\link{clusterGroupBound}()} to test all groups in a
  hierarchical clustering tree.
  Use \code{\link{groupBound}()} to compute the lower bound for selected
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

for (ii in unique(ind)) {
  id <- which(ind == ii)
  Sigma[id, id] <- rho
}
diag(Sigma) <- 1
print.table(Sigma, zero.print=".") ## depicting the 2 blocks

x <- matrix(rnorm(n * p), nrow = n) \%*\% chol(Sigma)

## Create response with active variables 1 and 21
beta    <- rep(0, p)
beta[1] <- 5

y  <- as.numeric(x \%*\% beta + rnorm(n))

## Compute the lower bound for all groups in a hierarchical clustering tree
cgb5 <- clusterGroupBound(x, y, nsplit = 5)

## Plot the tree with y-axis proportional to the (log) of the number of
## group members and node sizes proportional to the lower l1-norm bound.
plot(cgb5)

## Show the lower bound on the y-axis and node sizes proportional to
## number of group members
plot(cgb5, yaxis = "")
}
\keyword{htest}% from RShowDoc("KEYWORDS")
\keyword{regression}
