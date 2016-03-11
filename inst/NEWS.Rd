% Check from R:
% news(db= tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/hdi/pkg/inst/NEWS.Rd"))
\name{NEWS}
\title{hdi News}
\encoding{UTF-8}

\section{CHANGES IN VERSION 0.1-5)}{
  \itemize{
    \item new argument \code{suppress.grouptesting} in \code{ridge.proj}
  and \code{lasso.proj} to speed up computation.
    \item argument \code{verbose} instead of \code{trace} in
    \code{multi.split} and \code{stability}.
    \item new functions to generate reference data sets, see
    \code{rXb}, \code{rX} and \code{rb}.
    \item added NEWS and CITATION file.
  }
}
