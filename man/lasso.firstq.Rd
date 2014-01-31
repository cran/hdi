\name{lasso.firstq}
\alias{lasso.firstq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to determine the first q predictors in the lasso path.}
\description{Determines the q predictors that enter the lasso path first.}
\usage{lasso.firstq(x, y, q, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{q}{Number of predictors that should be selected.}
  \item{...}{Additional arguments to be passed to \code{\link{glmnet}}.}
}
\details{Function only calls glmnet in a special way, see source code.}
\value{Vector of selected predictors.}
\author{Lukas Meier}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{hdi}}.}
\examples{
x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)
sel <- lasso.firstq(x, y, q = 5)
sel
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}% __ONLY ONE__ keyword per line
