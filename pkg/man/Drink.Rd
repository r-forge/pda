\name{Drink}
\alias{Drink}
\title{Water Temperature Effect on Cow Drinking}
\usage{
data(Drink)
}
\format{
    Drink data frame with 261 observations on 6 variables.
    \tabular{rlll}{
	[,1] \tab cow \tab  factor \tab  cow identifier \cr
	[,2] \tab hc \tab   factor \tab  Heifer or Cow \cr
	[,3] \tab trt \tab  factor \tab  treatment group \cr
	[,4] \tab per \tab  factor \tab  period (Early/Middle/Late) \cr
	[,5] \tab dmi \tab  numeric \tab dry matter intake (DMI) \cr
	[,6] \tab coce \tab factor  \tab plot code
    }
}
\description{A dairy scientist is interested in the effect on milk yield of feeding
\code{cow}s hot (lukewarm, actually) instead of cold \code{water}.  This
may have economic importance if the temperature of \code{water} can
alter milk yield by even a pound per week.  Animals were put on hot
(or cold) \code{water} for three weeks, with measurements taken in the 
final week (as 7-day \code{milk} yield) of the \code{period}.  Each
\code{cow} was given both hot and cold \code{water} over a six week (two
\code{period}s), with \code{cow}s randomized as to whether they received
hot or cold \code{water} first in each pair.  \code{Cow}s might be
treated over several pairs of \code{period}s during the course of the
study. \code{Milk} yield should gradually decrease over time, regardless
of treatment.  This decline is confounded with the hot/cold treatment
for any given cow, but can be sorted out by comparing \code{cow}s given
hot or cold \code{first}. There is a covariate
\code{dim} (days in milk) that indicates how long the cow has
been producing milk; \code{milk} yield tends to rise initially
and then gradually fall, with a total lactation (milk producing)
time of roughly 305 days. In addition the \code{month} of entry
into the study is included to help assess seasonal effects if any.}
\source{Dave Combs}
\references{Wattiaux MA, Combs DK and Shaver RD (1994)
	`Lactational responses to ruminally undegradable protein by
	dairy cows fed diets based on alfalfa silage',
	\emph{J Dairy Science 77}, 1604-1617.}
\examples{
# Average measurement per cow
data( Drink )

Drink2 <- Drink[Drink$period < 3,]
Drink4 <- Drink[!is.na(match(Drink$cow,Drink$cow[Drink$period == 4])),]

Drink$cow <- factor( Drink$cow )
Drink$month <- ordered( Drink$month )
Drink$period <- ordered( Drink$period )
Drink.fit <- aov( milk ~ cow + period * water, Drink )

Drink2$cow <- factor( Drink2$cow )
Drink2$month <- ordered( Drink2$month )
Drink2$period <- ordered( Drink2$period )
Drink2.fit <- aov( milk ~ cow + period * water, Drink2 )

Drink4$cow <- factor( Drink4$cow )
Drink4$month <- ordered( Drink4$month )
Drink4$period <- ordered( Drink4$period )
Drink4.fit <- aov( milk ~ cow + period * water, Drink4 )

# I:27.1 Drink cross-over interaction plots
lsd.plot( Drink2.fit, factors = c("period","water"),
  xpos = 1.25, xlab = "(a) cows in first two periods",
  ylab = "7-day milk yield (lb)",
  main = "Figure 27.1",
  more = TRUE, split = c(1,1,2,1) )
lsd.plot( Drink4.fit, factors = c("period","water"),
  xpos = 1.75, xlab = "(b) cows in all four periods", ylab = "",
  main = "Drink",
  split = c(2,1,2,1) )

# better approach: mixed model fit using lme()
library(lme4)
Drink2.lme <- lmer(milk ~ period * water + (1|cow),
  data = Drink2 )
summary(Drink2.lme)
anova(Drink2.lme)
VarCorr( Drink2.lme)

int.plot( Drink2.lme, factors = c("period","water"),
   xpos = 1.25, xlab = "(a) cows in first two periods",
   ylab = "7-day milk yield (lb)", bar = "ellipse" )

Drink4.lme <- lmer(milk ~ period * water + (1|cow),
  data = Drink4 )
Drink4.lme
anova(Drink4.lme)
VarCorr( Drink4.lme)

int.plot( Drink4.lme, factors = c("period","water"),
   xpos = 1.75, xlab = "(b) cows in all four periods",
   ylab = "", bar = "ellipse" )
int.plot( Drink4.lme, factors = c("water","period"),
   xpos = 1.75, xlab = "(c) cold vs hot",
   ylab = "", bar = "ellipse" )
} 
% Converted by  version .
\keyword{datasets}
