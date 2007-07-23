\name{pvalue.ellipse}
\alias{pvalue.ellipse}
\title{Center and Radius for p-value Ellipses}
\description{
This is used by \code{int.plot} to construct ellipses for means at the
level of the x factor that are not significantly different.
Currently this does \emph{not} handle mixed models properly.
}
\usage{
pvalue.ellipse(plet, lsmeans, se, rdf,
   level=0.05, scale=0.3, full.width=TRUE, pset)
}
\arguments{
  \item{plet}{character vector with letter codes for non-significant differences}
  \item{lsmeans}{least squares means}
  \item{se}{standard errors for means}
  \item{rdf}{residual degrees of freedom}
  \item{level}{significance level}
  \item{scale}{scaling of ellipse around means }
  \item{full.width}{keep ellipses as full width to show how close means
    are if TRUE}
  \item{pset}{letters to use for means that are significantly different}
}
\details{
  The ellipse will be the height of an LSD comparison interval if
  \code{full.width} is TRUE, or otherwise will be \code{se*scale} beyond
  the outer means of a comparison group.
}
\value{
  A \code{2*n} matrix, where \code{n} is the number of distinct
  ellipses. The first row contains ellipse centers while the second row
  contains the radii.
}
\author{Brian S. Yandell}
\seealso{\code{\link{lsmean}}, \code{\link{pletter}},
  \code{\link{int.plot}}, \code{\link{ellipse}}}
  
\examples{
example( Oocyte )
int.plot( Oocyte.fit, bar.plot = "ellipse" )
}
\keyword{htest}
