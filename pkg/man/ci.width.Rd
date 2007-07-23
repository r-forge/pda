\name{ci.width}
\alias{ci.width}
\title{utility functions for ci.plot}
\usage{ci.width(object, data, factors, lsm,
  rdf=df.resid(object), level=0.05,  crit=qt(1-level/2, rdf), \dots)}
\arguments{
  \item{object}{fitted model object (optional unless response is
    omitted).}
  \item{data}{data frame in which to interpret the variables named
    in the fit \code{object} (taken from \code{object} if omitted).}
  \item{factors}{Character string of length 2 with names of
    \code{x.factor} and \code{trace.factor} as found in \code{object}
    and data (default is first 2 factors in \code{object}).}
  \item{lsm}{least squares means (default taken from object)}
  \item{rdf}{residual degrees of freedom}
  \item{level}{significance level for confidence interval}
  \item{crit}{value from t distribution at level with rdf degrees of
    freedom by default}
  \item{\dots}{ignored}
}
\description{}
\value{Width of confidence interval at all levels of factor combinations,
  computed as (critical.value*sd/size).}
\seealso{\code{\link{ci.plot}}.}
\keyword{htest}
% Converted by  version .
