\name{all.factors}
\alias{all.factors}
\alias{all.covars}
\alias{all.predictors}
\alias{model.response}
\title{determine names of all factors, covariates, predictors (factors \&
  covariates)or response used in model fit}
\usage{
all.factors(object, data)
all.covars(object,data)
all.predictors(object, data)
}
\arguments{
  \item{object}{object from linear model fit}
  \item{data}{data frame in which to interpret the variables named
    inthe full and reduced objects (determined from reduced fit if
    omitted)}}
\value{list of names of factors, covariates or predictors (factors \&
  covariates) or response used in model fit
}
\description{}
\seealso{\code{\link[base]{all.names}}.}
\examples{
\dontrun{
fit <- lm(y~x+a*b,data)
# return c("a","b")all.factors(fit)
# return "x"all.covars(fit)
# return c("x","a","b")all.predictors(fit)
# return "y"model.response(fit)
}
}
% Converted by  version .
\keyword{design}
