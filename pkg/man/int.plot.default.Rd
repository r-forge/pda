\name{int.plot.default}
\alias{ci.plot.default}
\alias{int.plot.default}
\title{Interaction Plots of Response with Interval Bars}
\usage{
\method{ci.plot}{default}(object, response, group, xlab, ylab,
  bar.plot="ci", \dots)
\method{int.plot}{default}(object, response, group, xlab, ylab, \dots)
}
\arguments{
  \item{object}{factor to be plotted on the horizontal axis.}
  \item{response}{numeric vector.}
  \item{group}{factor whose levels will be separate traces.}
  \item{\dots}{optional specification of other parameters to be passed
    to plotting routines and/or to the bar-drawing routines
    \code{lsd.bar} and \code{se.bar}. There are some parameters, such
    as \code{width} and \code{edge} which can change to plotting
    area.}
  \item{bar.plot}{"lsd" (default for \code{int.plot}), "ci" (default for
    \code{ci.plot}), "test", "ellipse" or "none". The "test" is a
    sqrt(2) contraction of confidence interval for approximate graphical
    comparisons of means.
  } 
  \item{xlab,ylab}{x and y axis labels}
}
\description{An interaction plot shows the mean responses for each
level of the \code{object} at each level of the \code{group}.
}
\value{list containing sd, rdf, level and sample size, and either lsd
    (lsd.plot) or width (ci.plot).}
  \seealso{\code{\link{int.plot}},
    \code{\link[base]{xyplot}}, \code{\link{margin.plot.lm}},
    \code{\link{lsd.bar}}, \code{\link{ci.width}},
    \code{\link{pvalue.ellipse}}.} 
\examples{
data(ToothGrowth)
attach(ToothGrowth)
lsd.plot(dose, len, supp)
detach()
}
\keyword{design}
% Converted by  version .
