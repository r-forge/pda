\name{int.plot.lm}
\alias{lsd.plot.lm}
\alias{ci.plot.lm}
\alias{ci.plot.lme}
\alias{ci.plot.lmer}
\alias{int.plot.lm}
\alias{int.plot.lme}
\alias{int.plot.lmer}
\title{Linear Model-driven Interaction Plots of Response with Interval Bars}
\usage{
\method{ci.plot}{lm}(object, \dots, bar.plot="ci")
\method{ci.plot}{lme}(object, data, \dots, bar.plot="ci")
\method{ci.plot}{lmer}(object, data, \dots, bar.plot="ci")
\method{int.plot}{lm}(object, data, factors, lsm, ci, offset,
   bar.plot="lsd", ylim, type, width, sort.mean=FALSE, edge,
   xlab, ylab, level, rdf, fine, xpos, ypos, cex, white, panelf = NULL,\dots)
\method{int.plot}{lme}(object, data, \dots)
\method{int.plot}{lmer}(object, data, \dots)
}
\arguments{
  \item{object}{fitted model object (optional unless response is
    omitted).}
  \item{data}{data frame in which to interpret the variables named
    in the fit \code{object} (taken from \code{object} if omitted).}
  \item{factors}{Character string of length 2 with names of
    \code{x.factor} and \code{trace.factor} as found in \code{object}
    and \code{data} (default is first 2 factors in \code{object}).}
  \item{lsm}{Least squares means for plotting (default is object from
    model).}
  \item{sort.mean}{should \code{x.factor} levels be sorted by mean value
    (default=FALSE).}
  \item{ylim,type,xlab,ylab}{plot parameters}
  \item{ci,offset,width,edge,xpos,ypos,level,rdf,fine,cex}{optional
    parameters for bar-drawing routines \code{lsd.bar} and
    \code{se.bar}. There are some parameters, such as \code{width} and
    \code{edge} which can change to plotting area.}
  \item{bar.plot}{"lsd" (default for \code{int.plot}), "ci" (default for
    \code{ci.plot}), "test", "ellipse" or "none". The "test" is a
    sqrt(2) contraction of confidence interval for approximate graphical
    comparisons of means.}
  \item{white}{white plot background if TRUE (default)}
  \item{panelf}{optional extra panel function for lattice plot}
  \item{\dots}{optional parameters for plotting routines}
}
\description{A plot is created showing the requested function of responses
  for each level of the \code{x.factor} at each level of the
  \code{trace.factor}.
}
\value{list containing sd, rdf, level and sample size, and either lsd
    (lsd.plot) or width (ci.plot).}
  \seealso{\code{int.plot},
    \code{\link[lattice]{xyplot}}, \code{\link{margin.plot.lm}},
    \code{\link{lsd.bar}}, \code{\link{ci.width}},
    \code{\link{pvalue.ellipse}}.} 
\examples{
data(ToothGrowth)
ToothGrowth$dose = ordered(ToothGrowth$dose)
tooth.fit = aov(len~dose*supp, ToothGrowth)
lsd.plot( tooth.fit, ToothGrowth, factors = c("supp","dose"))

cat.lsm = lsmean(tooth.fit)
lsd.plot(tooth.fit, ToothGrowth, factors = c("supp","dose"),
         lsm = cat.lsm)
}
\keyword{design}
% Converted by  version .
