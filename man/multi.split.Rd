\name{multi.split}
\alias{multi.split}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to calculate p-values based on multi-splitting approach}

\description{Function to calculate p-values based on multi-splitting approach}
\usage{
multi.split(x, y, B = 100, fraction = 0.5, ci = TRUE, ci.level = 0.95,
            model.selector = lasso.cv,
            classical.fit = lm.pval, classical.ci = lm.ci,
            parallel = FALSE, ncores = 4,
            gamma = seq(ceiling(0.05*B)/B,1-1/B,by=1/B),
            args.model.selector = NULL, args.classical.fit = NULL,
            args.classical.ci = NULL,
            return.nonaggr = FALSE, return.selmodels = FALSE, trace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Design matrix (without intercept).}
  \item{y}{Response vector.}
  \item{B}{Number of sample-splits.}
  \item{fraction}{Fraction of data used at each sample split for
    the model selection process. The remaining data is used for
    calculating the p-values.}
  \item{ci}{Should a confidence interval be calculated for every
    parameter? (logical)}
  \item{ci.level}{Coverage level of confidence interval.}
  \item{model.selector}{Function to perform model selection. Default is
    \code{\link{lasso.cv}}. Function must have at
    least two arguments: x (the design matrix) and y (the response
    vector). Return value is the index vector of selected columns. See
    \code{\link{lasso.cv}} and \code{\link{lasso.firstq}} for an 
    example. Additional arguments can be passed through
    \code{args.model.selector}.}
  \item{classical.fit}{Function to calculate (classical)
    p-values. Default is \code{\link{lm.pval}}. Function must 
    have at least two arguments: x (the design matrix) and y (the
    response vector). Return value is the vector of p-values. See
    \code{\link{lm.pval}} for an example. Additional arguments can be
    passed through \code{args.classical.fit}.}
  \item{classical.ci}{Function to calculate (classical)
    confidence intervals. Default is \code{\link{lm.ci}}. Function must 
    have at least 3 arguments: x (the design matrix), y (the
    response vector) and level (the coverage level). Return value is the
    matrix of confidence intervals. See \code{\link{lm.ci}} for an
    example. Additional arguments can be passed through
    \code{args.classical.ci}.}
  \item{parallel}{Should parallelization be used? (logical)}
  \item{ncores}{Number of cores used for parallelization.}
  \item{gamma}{Vector of gamma-values. In case gamma is a scalar, the
    value Q_j instead of P_j is being calculated (see reference below).}
  \item{args.model.selector}{Named list of further arguments for
    function \code{model.selector}.}
  \item{args.classical.fit}{Named list of further arguments for function
    \code{classical.fit}.}
  \item{args.classical.ci}{Named list of further arguments for function
    \code{classical.ci}.}
  \item{return.nonaggr}{Should the unadjusted p-values be reported? (logical).}
  \item{return.selmodels}{Should the selected models (at each split) be
    reported? (logical).}
  \item{trace}{Should information be printed out while computing? (logical).}
}
\details{...}
\value{
  \item{pval.corr}{Vector of multiple testing corrected p-values.}
  \item{gamma.min}{Value of gamma where minimal p-values was attained.}
}
\references{
  Meinshausen, N., Meier, L. and B\"uhlmann, P. (2009),
  \emph{P-values for high-dimensional regression}, Journal of the
  American Statistical Association 104, 1671-1681.
}
\author{Lukas Meier}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{...}

\examples{
x <- matrix(rnorm(100*1000), nrow = 100, ncol = 1000)
y <- x[,1] * 2 + x[,2] * 2.5 + rnorm(100)

## Multi-splitting with lasso.firstq as model selector function
fit.multi <- multi.split(x, y, model.selector =lasso.firstq,
                         args.model.selector = list(q = 10))
fit.multi
fit.multi$pval.corr[1:10] ## the first 10 p-values
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression} % __ONLY ONE__ keyword per line
