\name{sample.size}
\alias{sample.size}
\alias{sample.size.default}
\alias{sample.size.lm}
\title{compute sample size by group}
\usage{
sample.size(object, \dots)
\method{sample.size}{default}(object, \dots)
\method{sample.size}{lm}(object, data, factors, na.action, \dots)
}
\arguments{
  \item{object}{vector or list of vectors, or object from linear model fit}
  \item{\dots}{additional vectors in implicit list}
  \item{data}{data frame in which to interpret the variables named
    inthe full and reduced objects (determined from reduced fit if
    omitted)}
  \item{factors}{Character string of length 2 with names of x.factor
    and trace.factoras found in fit and data (default is first 2
    factors in fit).}
  \item{na.action}{missing data action}
}
\value{In the case of a list of vectors, a vector of their lengths. For
  a linear model object, the sample size by factors.
}
\description{}
\examples{
\dontrun{
sample.size(list(x1,x2,x3))
sample.size(x1,x2,x3)
sample.size(lm(y~a+b,data))
}
}
\keyword{htest}
% Converted by  version .
