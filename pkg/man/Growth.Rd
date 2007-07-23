\name{Growth}
\alias{Growth}
\title{Growth Medium for Plant Tissue}
\usage{
data(Growth)
}
\format{
    Growth data frame with 183 observations on 5 variables.
    \tabular{rlll}{
	[,1] \tab trt \tab  factor \tab  treatment identifier \cr
	[,2] \tab BA \tab  factor \tab  BA level \cr
	[,3] \tab TDZ \tab  factor \tab  TDZ level \cr
	[,4] \tab code \tab  factor \tab plot code \cr
	[,5] \tab advplt \tab  numeric \tab number of adventitious plants
    }
}
\description{This research concerned growth medium for cloning plants
(Stieve, Stimart and Yandell 1992).  Plant tissue (a small part of a
leaf) was placed in a Petri dish which had a `food preparation' on it.
The dish was sealed and placed in a controlled chamber with adequate
light and proper temperature.  The `food preparation', or medium, was
designed to encourage growth of new plant shoots from the tissue.  The
response of interest here was the number of `adventitious shoots'
(\code{advplt}), or new plants growing out of the tissue.  Basically,
the more adventitious shoots produced, the better the growing media.
The adventitious plants can be separated and grown into numerous
individuals for greenhouse cultivation, saving tremendous space and
cost for nurseries and florists. The key question is: what is the best
growing condition for the shoots?

There were 12 factor combinations (\code{trt}) of interest,
coming from a two-factor layout with added control (BHTA, \code{trt} =
20).  One factor combination was missing; the plant material seemed
to be burned by high levels of both chemicals together.
}
\source{Sue Stieve \& Dennis Stimart}
\references{Stieve SM, Stimart DP and Yandell BS (1992)
        `Heritable tissue culture induced variation in \emph{Zinnia
        Marylandica}',
        \emph{Euphytica 64}, 81-89.}
\seealso{\code{\link{GrowthCount}}.}
\examples{
data( Growth )

#Growth <- Growth[,c("trt","advplt")]
# drop the added BHTA control ( trt = 20 )
Growth <- Growth[ Growth$trt != 20, ]
Growth$logadvplt <- log10(1+Growth$advplt)
 
Growth$trt <- factor(Growth$trt)
Growth$BA <- ordered(Growth$BA)
Growth$TDZ <- ordered(Growth$TDZ)
Growth <- Growth[!is.na(Growth$advplt),]

sample.size( Growth$BA, Growth$TDZ )

Growth.zero <- Growth[Growth$advplt>0,]

Growth.fit <- lm(logadvplt~trt,Growth)
Growth.zero.fit <- lm(logadvplt~trt,Growth.zero)

Growth.full <- aov( logadvplt ~ BA * TDZ, Growth )
Growth.add <- aov( logadvplt ~ BA + TDZ, Growth )

# D:11.1 Growth interaction plots with a missing cell

int.plot( Growth.full, xpos = 2.5, ypos = .425,
  xaxt = "n", xlab = "(a) BA by TDZ", ylab = "",
  main = "Figure D:11.1",
  more = TRUE, split = c(1,1,2,1) )
int.plot( Growth.full, fac = c("TDZ","BA"),
  xpos = 3.25, ypos = .425, xaxt = "n",
  xlab = "log(1 + advplt)", ylab = "(b) TDZ by BA",
  main = "Growth Interaction",
  split = c(2,1,2,1) )

# E:13.2 Growth residual plots with zeros
tmpdata = data.frame( x = jitter( predict( Growth.fit )),
  y = jitter( resid( Growth.fit )), group = Growth$code )
print( xyplot( y ~ x, tmpdata, groups = group,
  panel = function(...) {
    panel.superpose(...)
    panel.abline( h = 0, lty = 2 )
    panel.abline( h = c(-1,1) * std.dev( Growth.fit ), lty = 3 )
    panel.abline( 0, -1 )
  },
  xlab = "(a) mean with zeroes", ylab = "residuals",
  main = "Figure E:13.2" ),
  more = TRUE, split = c(1,1,2,1) )

tmpdata = data.frame( x = jitter( predict( Growth.zero.fit )),
  y = jitter( resid( Growth.zero.fit )), group = Growth.zero$code )
print( xyplot( y ~ x, tmpdata, groups = group,
  panel = function(...) {
    panel.superpose(...)
    panel.abline( h = 0, lty = 2 )
    panel.abline( h = c(-1,1) * std.dev( Growth.zero.fit ), lty = 3 )
    panel.abline( 0, -1 )
  },
  xlab = "(b) mean without zeroes", ylab = "",
  main = "Growth Residuals" ), 
  split = c(2,1,2,1) )
} 
% Converted by  version .
\keyword{datasets}
