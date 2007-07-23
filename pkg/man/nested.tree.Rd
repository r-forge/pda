\name{nested.tree}
\alias{nested.tree}
\title{draw nested tree}
\description{draw nested tree based on model estimates/predictors}
\usage{
nested.tree(object, nesting=names(means),
means=nested.means(object, nesting[-1]), offset=0.05, small=0.02,
estimate="estimate", bar="stderr", jitter=offset/2, \dots)}
\arguments{
    \item{object}{data frame}
    \item{nesting}{character vector of nesting names for \code{means}}
    \item{means}{list of means with names given by \code{nesting}}
    \item{offset}{offset for accuracy \code{bar}}
    \item{small}{smallest relative separation of names allowable
	(scalar, or vector withsame length as \code{nesting})}
    \item{estimate}{name of estimate or prediction column}
    \item{bar}{name of column for accuracy}
    \item{jitter}{jitter of names}
    \item{\dots}{plotting parameters}}
\value{}
\details{Either (\code{nesting} and \code{object}) or \code{means} must be
    specified.}
\examples{
\dontrun{
nested.tree(means=wheatfix.means,bar="ci",small=c(0,0,.05))
}
}
% Converted by  version .
\keyword{hplot}
