\name{nested}
\alias{nested}
\alias{nested.default}
\alias{nested.aovprojlist}
\title{nested weighted means for experimental units}
\description{This is under development.  It is a simple approximation to
    reducemeasurements to a single value for each sized experimental
    unit.}
  \usage{
nested( object, \dots)
\method{nested}{default}(object, data, project = proj(object), \dots)
\method{nested}{aovprojlist}( object, data, factors, response="y",
  nesting=c(1,ncol(ref)), ref = projref( object ), \dots)
}
\arguments{
    \item{object}{linear model fit object with Error() nesting or
      projection list}
    \item{data}{data frame used in object}\item{project}{projection of data
	using object}
    \item{factors}{vector of names of factors in
	object}\item{response}{name of response}
    \item{nesting}{levels of nesting to be averaged (1=intercept,
	2=whole plot, 3=sub-plot,etc.)}
    \item{ref}{data frame of reference values from
      projection}
    \item{\dots}{additional arguments for \code{projwt}}}
\value{Data frame with average \code{response} for unique \code{factor}
    combinations and additional column \code{weight} with counts of
    cases used in averages.}
\details{See \code{nested.default} and \code{nested.aovprojlist} for
    code details.}
\seealso{\code{\link[base]{proj}}.}
\examples{
\dontrun{
nest.data <- nested( proj(object), data, c("A","b"), "y", c(1,2))
aov(y~A*B, nest.data)
}}
% Converted by  version .
\keyword{design}
