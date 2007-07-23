\name{Hardy}
\alias{Hardy}
\title{Cold Hardiness Adaptation}
\usage{
data(Hardy)
}
\format{
    Hardy data frame with 75 observations on 6 variables.
    \tabular{rlll}{
	[,1] \tab potato \tab  factor \tab potato identifier \cr
	[,2] \tab regime \tab  factor \tab hardiness regime \cr
	[,3] \tab temp \tab  factor \tab  temperature \cr
	[,4] \tab code \tab  factor \tab plot code \cr
	[,5] \tab photo \tab  numeric \tab photosynthesis score \cr
	[,6] \tab leak \tab  numeric \tab ion leakage score
    }
}
\description{Scientists in horticulture have been interested in how plants adapt
to cold climates such as we find in Wisconsin.  They have noticed that
plants which have been previously `acclimatized' suffer less damage.
Palta and Weiss (1993) conducted a
controlled laboratory experiment to examine two species of
\code{potato}es (1, 2) subjecting plants to one of two acclimatization
\code{regime}s (1 kept in cold room; 0 kept at room temperature) for
several days.  Plants were later subjected to one of two cold
\code{temp}erature treatments (1 just below freezing, \emph{-4 C}; 2 way below
freezing, \emph{-8 C}).  Two responses were measured, both being damage
scores for \code{photo}synthesis and ion \code{leak}age.
Initially there were 80 plants, but some were lost during the experiment.
Your analysis should take care to address the imbalance due to unequal
sample sizes for the different factor combinations.
}
\source{Jiwan P Palta (\url{mailto:palta@calshp.cals.wisc.edu})
  \& Laurie S Weiss, U WI Horticulture, (\url{http://www.hort.wisc.edu})
}
\references{Palta JP and Weiss LS (1993)
	`Ice formation and freezing injury: an overview on the
	survival mechanisms and molecular aspects of injury and cold
	acclimation in herbaceous plants',
	in \emph{Advances in Plant Cold Hardiness},
	ed. by PH Li and L Christersson.
	CRC Press, Boca Raton, LA.}
\examples{
data( Hardy )

Hardy$potato <- factor(Hardy$potato)
Hardy$regime <- factor(Hardy$regime)
Hardy$temp <- factor(Hardy$temp)

# full model: you may want to consider various reduced models */
Hardy.fit <- aov(leak ~ potato * regime * temp, Hardy)
# Type I sums of squares
summary( Hardy.fit )
# Type III sums of squares
drop1( Hardy.fit, formula( Hardy.fit ), test = "F" )

# Least Squares means
Hardy.lsm <- lsmean(Hardy.fit)
# default plots for lm object
plot(Hardy.fit)

# diagnostic plot of resid vs predict with 0 line and +/- 1 SD
Hardy$code = factor( Hardy$code )
xyplot( resid(Hardy.fit) ~ jitter(predict(Hardy.fit)),
  groups = Hardy$code, pch = levels( Hardy$code ), cex = 1.5,
  panel = function(x,y,sd=std.dev(Hardy.fit),...){
    panel.xyplot(x,y,...)
    panel.abline(0,0,lty=2)
    panel.abline(sd,0,lty=3)
    panel.abline(-sd,0,lty=3)
  } )

# interaction plot with LSD bars (but play with this!)
lsd.plot(Hardy.fit,Hardy,c("potato","regime"),
  more = TRUE, split = c(1,1,2,1) )
# the following give a slightly different plot -- why?
attach(Hardy)
lsd.plot(potato,leak,regime, split = c(2,1,2,1))
detach()

# side-by-side interaction plots for the subsets of t=1,2
lsd.plot(Hardy.fit,Hardy,c("potato","regime"),
  Hardy.lsm[Hardy.lsm$temp==1,],
  more = TRUE, split = c(1,1,2,1) )
lsd.plot(Hardy.fit,Hardy,c("potato","regime"),
  Hardy.lsm[Hardy.lsm$temp==2,],
  split = c(2,1,2,1) )
} 
% Converted by  version .
\keyword{datasets}
