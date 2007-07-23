\name{df.resid}
\alias{df.resid}
\alias{df.resid.default}
\alias{df.resid.lm}
\alias{df.resid.lme}
\alias{df.resid.lmer}
\title{residual degrees of freedom}
\usage{
df.resid(object, \dots)
\method{df.resid}{default}(object, \dots)
\method{df.resid}{lm}(object, resid=residuals(object), \dots)
\method{df.resid}{lme}(object, resid = residuals( object ), factors, \dots)
\method{df.resid}{lmer}(object, \dots)
}
\arguments{
  \item{object}{vector, list of vectors, or object from linear model
    fit}
  \item{resid}{Residuals of \code{object}}
  \item{factors}{Character string of factors for \code{object}}
  \item{\dots}{more vectors}
}
\description{}
\value{If a single vector is given, return its degrees of freedom
  corrected for the mean, length(object)-1. If a list of vectors is
  given (explicitly or implicitly), sum up df for each. If object is a
  linear model fit, extract its \code{df.resid} element or make
  appropriate guess based on model degrees of freedom if \code{df.resid}
  is missing. The \code{lmer} method attempts to find total random
  effects degrees of freedom.}
\examples{
\dontrun{
df.resid(x)df.resid(lm(y~a+b,data))
}
}
\keyword{htest}
% Converted by  version .
