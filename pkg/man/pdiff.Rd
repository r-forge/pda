\name{pdiff}
\alias{pdiff}
\title{Table of p-values for Differences of Least Squares Means}
\description{
Takes the fit of a linear model, extracts the least squares means and
their standard errors, and the residual DF, and constructs the p-values
for all pairwise comparisons.
}
\usage{
pdiff(fit, lsm=lsmean(fit), pred, se, rdf=df.resid(fit))
}
\arguments{
 \item{fit}{object from linear model fit}
 \item{lsm}{least squares means object from \code{lsmean}}
 \item{pred}{outer difference of least squares means}
 \item{se}{SE for differences in \code{pred}}
 \item{rdf}{residual degrees of freedom from model fit}
}
\value{
Square matrix of size \code{length( pred )} sorted from largest to
smallest \code{pred} with p-values as
entries. Diagonal entries are all 1. Row and column names are levels of \code{pred}.
}
\author{Brian S. Yandell}
\seealso{\code{\link{lsmean}}, \code{\link{pletter}}, \code{\link{pvalue.ellipse}}}

\examples{
example( Oocyte )
pdiff( Oocyte.fit )
}
\keyword{htest}
