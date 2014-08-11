\name{clusterGroupBound}
\alias{clusterGroupBound}
\title{
 Group test of variable importance in a high-dimensional linear model,
 using a hierarchical structure.} 
\description{
 Computes confidence intervals for the l1-norm of groups of regression
 parameters in a hierarchical clustering tree.}
\usage{
  clusterGroupBound(x, y, method = "average",
                    dist = as.dist(1 - abs(cor(x))), alpha = 0.05,
                    hcloutput, nsplit = 11, s = min(10, ncol(x) - 1),
                    silent = FALSE, setseed = TRUE, lpSolve = TRUE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 The design matrix of the regression with p columns for p predictor
 variables and n rows that correspond to n observations.}
  \item{y}{
 The response variable; a numeric vector of length n.}
  \item{method}{
 The method used for constructing the hierarchical clustering tree
 (default is "average" linkage). Alternatively, you can provide your own
 hierarchical clustering through the optional argument \code{hcloutput}.}
  \item{dist}{
 A distance matrix can be entered as an argument, on which the
 hierarchical clustering will be based. The default option is that the distance
 between variables will be calculated as 1 less the absolute correlation
 matrix. Alternatively, you can provide your own hierarchical clustering
 through the optional argument \code{hcloutput}.}
  \item{alpha}{
 The level in (0, 1) at which the confidence intervals are to be
 constructed.}
\item{hcloutput}{
  Optional argument.
  The output of a call the the hclust function. If it is provided, the
  arguments dist and method are ignored.}
 \item{nsplit}{
 The number of data splits used.}
  \item{s}{
 The dimensionality of the projection that is used. Lower values lead
 to faster computation and if n > 50, then s is set to 50 if left
 unspecified to avoid lengthy computations.}
\item{silent}{
 Output is supressed if this option is set to true.}
  \item{setseed}{
 If \code{setseed} is true (recommended), then the same random seeds are
 used for all groups, which makes the confidence intervals
 simultaneously valid over all groups of variables tested.}
  \item{lpSolve}{
 Only set to false if lpSolve is not working on the current machine
 (setting it to false will results in much slower computations; only use
 on small problems).}
}
%\details{
%}
\value{
Returns a list with components
 \item{groupNumber}{The index of the group tested in the original
   hierarchical clustering tree} 
 \item{members}{A list containing the variables that belong into each
   testes group} 
\item{noMembers}{A vector containing the number of members in each group}
\item{lowerBound}{The lower bound on the l1-norm in each group}
\item{position}{The position on the x-axis of each group (used for plotting)}
\item{leftChild}{Gives the index of the group that corresponds to the
  left child node in the tested tree (negative values correspond to leaf
  nodes)} 
\item{rightChild}{Same as \code{leftCHild} for the right child of each node}
\item{isLeaf}{Logical vector. Is \code{TRUE} for a group if it is a leaf
  node in the tested tree or if both child nodes have a zero lower bound
  on their group l1-norm} 
}
\references{
 Nicolai Meinshausen (2013) 
 Assumption-free confidence intervals for groups of variables in sparse
 high-dimensional regression. 
 http://arxiv.org/abs/1309.3489
}
\author{
 Nicolai Meinshausen meinshausen@stat.math.ethz.ch
}
%\note{
%}

\seealso{
  Use \code{clusterGroupBound} to test all groups in a hierarchical
  clustering tree.
  Use \code{groupBound} to compute the lower bound for selected
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

out <- clusterGroupBound(x, y, nsplit = 5)

## Plot and print the hierarchical group-test
plot(out)
print(out)
}
\keyword{confidence intervals}
\keyword{regression}
\keyword{hierarchical clustering}

