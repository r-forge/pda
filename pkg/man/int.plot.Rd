\name{int.plot}
\alias{lsd.plot}
\alias{ci.plot}
\alias{int.plot}
\title{Generic Interaction Plots of Response with Interval Bars}
\usage{
lsd.plot(object, \dots)
ci.plot(object, \dots)
int.plot(object, \dots)
}
\arguments{
  \item{object}{object to be plotted}
  \item{\dots}{optional specification of other parameters to be passed
    to methods, plotting routines and/or to the bar-drawing routines
    \code{lsd.bar} and \code{se.bar}. There are some parameters, such
    as \code{width} and \code{edge} which can change plotting
    area.
  }
}
\description{The routines \code{lsd.plot} simply calls
  \code{int.plot}, which is a method. \code{ci.plot}  is a method. An
  interaction plot is created. By 
default, lines for each value of the \code{group} are drawn in
differentline styles ('lsd.plot' and 'int.plot') and with plot symbols
so that they may be more easily distinguished.  This corresponds to plot
\code{type="b"}, which can be overridden. RA Fisher's LSD bar
(\code{lsd.plot}), confidence intervals (\code{ci.plot}) or test
ellipses may be added with the \code{bar.plot} option.
}
\value{list containing sd, rdf, level and sample size, and either lsd
    (lsd.plot) or width (ci.plot).}
\seealso{\code{\link[lattice]{xyplot}}, \code{\link{margin.plot}},
    \code{\link{lsd.bar}}, \code{\link{ci.width}}, \code{\link{int.plot.lm}},
    \code{\link{pvalue.ellipse}}.} 
\examples{
# the first use mimics interaction.plot
data(ToothGrowth)
ToothGrowth$dose = ordered(ToothGrowth$dose)
attach(ToothGrowth)
lsd.plot(dose, len, supp)
detach()
# or you could get more fancy by using a model fit
# here is an interaction plot based on linear model fit
ToothGrowth$dose = ordered(ToothGrowth$dose)
tooth.fit = aov(len~dose*supp, ToothGrowth)
lsd.plot( tooth.fit, ToothGrowth, factors = c("supp","dose"))
}
\keyword{design}
% Converted by  version .
