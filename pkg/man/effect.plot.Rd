\name{effect.fit}
\alias{effect.fit}
\alias{effect.fit.default}
\alias{effect.fit.listof}
\alias{effect.plot}
\title{Effects Plot}
\description{Get fit effects and plot beside residuals.
  Visually assess MS as spread of effects. See Emerson chapter in
  Hoaglin, Mosteller and Tukey (1991, ch. 8C) Fundamentals of
  Exploratory Analysis of Variance (Wiley).}
\usage{
effect.fit(object, \dots)
\method{effect.fit}{default}(object, data, factors, adjust=TRUE, influence=FALSE,
  terms.fit, \dots)
\method{effect.fit}{listof}(object, data, factors, terms.fit, stratum, \dots)
effect.plot(object, effect, adjust=TRUE, boxplot.min=10,
  jitter=0.1, width=jitter, xlim, ylim, xaxt, xlab, ylab, \dots)
}
\arguments{
  \item{object}{linear model (lm) fit}
  \item{data}{Data used in \code{object} fit}
  \item{factors}{Character string of factor names}
  \item{influence}{Include influence if \code{TRUE}}
  \item{terms.fit}{Character string of terms in \code{object} fit}
  \item{stratum}{Identifier of stratum to examine}
  \item{effect}{effect fit of \code{object}}
  \item{adjust}{adjust effects so that the average of squared effects
    is the mean squarefor that term}
  \item{boxplot.min}{use boxplot instead if more than this many values}
  \item{jitter}{uniform sideways jitter of plotted points
    (default=.1)}
  \item{width}{width of boxplot (default=jitter)}
  \item{xlim}{Plot parameters}
  \item{ylim}{Plot parameters}
  \item{xaxt}{Plot parameters}
  \item{xlab}{Plot parameters}
  \item{ylab}{Plot parameters}
  \item{\dots}{Plot parameters}
}
\value{list of effects}
\seealso{\code{\link{effect.plot}}.}
\examples{
data(ToothGrowth)
ToothGrowth$dose = ordered(ToothGrowth$dose)
effect.plot(lm(len~dose*supp, ToothGrowth))
}
\keyword{htest}
% Converted by  version .
