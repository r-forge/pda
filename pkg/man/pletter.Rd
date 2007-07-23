\name{pletter}
\alias{pletter}
\title{Letters for Means that are not Significantly Different}
\description{
Uses pdiff table to construct letters for multiple comparisons among means.
}
\usage{
pletter(pdiff, level=0.05, let=letters)
}
%- maybe also `usage' for other objects documented here.
\arguments{
 \item{pdiff}{table of p-values as from \code{pdiff()}}
 \item{level}{significance level }
 \item{let}{letters to use for non-significant differences}
}
\value{
Character vector with same length as \code{nrow( pdiff )} containing the
letter codes for the ordered means.
}
\author{Brian S. Yandell}
\seealso{\code{\link{lsmean}}, \code{\link{pdiff}}, \code{\link{pvalue.ellipse}}}

\examples{
example( Oocyte )
pletter( pdiff( Oocyte.fit ))
}
\keyword{htest}


