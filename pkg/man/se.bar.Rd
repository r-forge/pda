\name{se.bar}
\alias{se.bar}
\alias{lsd.bar}
\alias{lsd.bar.default}
\alias{lsd.bar.lm}
\alias{lsd.bar.lme}
\alias{lsd.bar.lmer}
\title{Standard Error and Fisher's LSD Bars for Interaction Plots}
\usage{
se.bar(xpos, ypos, bar, mod = "", cap = "SE", adj = 0, width = .01,
   lty, horiz, \dots)
lsd.bar(object, \dots)
\method{lsd.bar}{default}(object, response, group, \dots)
\method{lsd.bar}{lm}(object, data, factors, lsm, xpos, ypos,
  rdf = df.resid(object), level = .05, crit,
  cap = paste(100*level,"\%\nLSD",sep=""), mod, tol, \dots)
\method{lsd.bar}{lme}(object, data, \dots)
\method{lsd.bar}{lmer}(object, data, \dots)
}
\arguments{
  \item{xpos}{x position(s) of bar(s).}
  \item{ypos}{y position(s) of bar(s).}
  \item{bar}{height(s) of bar(s).}
  \item{mod}{modifier to bar caption ("*" used in lsd.bar for unequal
    sample sizes)}
  \item{cap}{caption for each bar (default for LSD with level; sets up standard
    deviation (SD) if set to "SD"; omitted if "")}
  \item{adj}{0 = left adjust caption to right of bar,1 = right adjust
    caption to left of bar.}
  \item{width}{relative width of bar ends (can be 0).}
  \item{object}{factor to be plotted on the horizontal axis (for
    \code{lsd.bar.default}) or fitted model object (for
    \code{lsd.bar.lm}).}
  \item{lsm}{least squares means object from \code{lsmean}}
  \item{group}{factor whose levels will be separate traces.}
  \item{response}{numeric vector.}
  \item{data}{data.frame in which to interpret the variables named
    in \code{object} (taken from object if omitted).}
  \item{factors}{Character string of length 2 with names of \code{object}
    and \code{trace.factor} as found in object and data (default is first 2
    factors in \code{object}).}
  \item{rdf}{residual degrees of freedom.}
  \item{level}{significance level.}
  \item{crit}{2-sided tail critical value of t-distribution.}
  \item{tol}{tolerance level}
  \item{lty}{line types}
  \item{horiz}{put bar horizontal if TRUE}
  \item{\dots}{optional specification of other parameters to be passed
    to se.bar. Standard error bars (se.bar) or an LSD bar (lsd.bar)
    are plotted at position(s) \code{cbind(xpos,ypos)}.  Lsd.bar uses
    the S method, with the \code{lsd.bar.default} calling
    \code{lsd.bar.lm} after setting up the fit, data and factor.
    \code{lsd.bar.lm} in turn calls \code{se.bar}
    after determining the barheight using sample sizes, rdf, sd and level.
    Confidence intervals or othermultiple bars can be constructed
    \code{ci.bar} (see \code{ci.plot}).}}
\value{a list containg xpos, ypos, bar, mod, cap, adj, width.}
\description{
Use \code{se.bar} to see the standard error. Use \code{lsd.bar} inside a
Lattice panel to display an LSD bar. Setting \code{cap="SD"} displays an
SD rather than LSD. Setting \code{crit=1} displays an SE. Asterisk
\code{*} beside bar indicated that harmonic mean of unequal sample sizes
was used in computation.}
\seealso{\code{\link{lsd.plot}}, \code{\link{ci.plot}}.}
\examples{
data(ToothGrowth)
attach(ToothGrowth)
lsd.bar(dose, len, supp)
detach()
}
\keyword{design}
% Converted by  version .
