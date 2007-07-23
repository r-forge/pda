\name{factor.labels}
\alias{factor.labels}
\title{determine labels for all factors used in model fit}
\usage{
factor.labels(object, data)
}
\description{}
\arguments{
  \item{object}{object from linear model fit}
  \item{data}{data frame in which to interpret the variables named
    inthe full and reduced objects (determined from reduced fit if
    omitted)}
}
\value{list containing labels of factors used in model
  fit}
\seealso{\code{\link{lsd.plot}},\code{\link{margin.plot.lm}}.}
\examples{
\dontrun{
# return list of length 2 with labels for a and b
factor.labels(lm(y~a+b,data))
}
}
% Converted by  version .
\keyword{design}
