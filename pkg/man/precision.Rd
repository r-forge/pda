\name{precision}
\alias{precision}
\title{round value to precision}
\usage{precision(x,se,digits=0)}
\description{}
\arguments{
  \item{x}{numeric vector}
  \item{se}{standard error(s) of value(s).}
  \item{digits}{one less than the number of extra digits}
  }
\value{value(s) rounded off to precision of 2 standard errors plus
  (\code{digits}+1) extra digit.}
\examples{
precision(75.111,.12)
}
\keyword{}
% Converted by  version .
