\name{Forage}
\alias{Forage}
\alias{ForCow}
\title{Forage Effect on Dry Matter Intake}
\usage{
data(ForCow)
data(Forage)
}
\format{
    ForCow data frame with 87 observations on 5 variables.
    \tabular{rlll}{
	[,1] \tab hc \tab  factor \tab  Heifer or Cow \cr
	[,2] \tab trt \tab  factor \tab  treatment group \cr
	[,3] \tab code \tab factor  \tab plot code \cr
	[,4] \tab id \tab  factor \tab  cow identifier \cr
	[,5] \tab dmi \tab  numeric \tab dry matter intake (DMI)
    }
    Forage data frame with 261 observations on 6 variables.
    \tabular{rlll}{
	[,1] \tab cow \tab  factor \tab  cow identifier \cr
	[,2] \tab hc \tab   factor \tab  Heifer or Cow \cr
	[,3] \tab trt \tab  factor \tab  treatment group \cr
	[,4] \tab per \tab  factor \tab  period (Early/Middle/Late) \cr
	[,5] \tab dmi \tab  numeric \tab dry matter intake (DMI) \cr
	[,6] \tab coce \tab factor  \tab plot code
    }
}
\description{The data for this problem are abstracted from a research project at
the US Center for Dairy Forage.  The experiment involved comparing the
effects of five different \code{diet}s on dry matter intake
(\code{dmi}), the amount of food eaten.  The animals under study were
heifers (cows with their first offspring) and mature cows during the
lactation (milk-producing) cycle following birth of a calf.  The five
diets ranged from low (\code{diet=1}) to high (\code{diet=5}) alfalfa
content. The experiment was actually fairly complicated.  The proportion
of alfalfa for the low group began at 45\% and increased to 65\% over
course of the lactaction cycle, while the high began at 85\% and
increased to 95\%.  However, the order of the groups (low to high)
remained the same. Interest here focuses on dry matter intake for
cows.}
\source{Tilak Dhiman}
\references{Dhiman TR, Kleinmans J, Tessmann NJ, Radloff HD and LD Satter (1995)
	Digestion and energy balance in lactating dairy cows fed
	varying ratios of alfalfa silage and grain,
	\emph{J Dairy Science 78}, 330.

	Dhiman TR and Satter LD (1996)
	Rumen degradable protein and its effect on microbial protein
	synthesis,
	\emph{J Dairy Science 00}.}
\examples{
# Average measurement per cow
data( ForCow )

# trt is ordered factor
ForCow$trt <- ordered( ForCow$trt )

# I do heifer, you do cow
ForCowh <- ForCow[ForCow$hc=="H", ]

# stem-and-leaf plot of raw data
stem( ForCowh$dmi, 5 )
# box-plots by treatment
print( bwplot( dmi ~ trt, ForCowh ),
  main = "Forage Raw Data by trt" )

# analysis of variance for heifer
ForCowh.aov <- aov( dmi ~ trt, ForCowh )
ForCowh.aov
summary( ForCowh.aov )

# linear model for heifer -- same idea
ForCowh.lm <- lm( dmi ~ trt, ForCowh )
ForCowh.lm
summary( ForCowh.lm )		# linear models summary
summary.aov( ForCowh.lm )		# anova-type summary

# ordinary means
tapply( ForCowh$dmi, ForCowh$trt, mean )
# predicted values for all responses
data.frame( trt=ForCowh$trt, mean=predict( ForCowh.aov ) )
# predicted values and their standard errors by trt
ForCowh.mean <- predict( ForCowh.aov, se=TRUE )
ForCowh.mean <- data.frame( trt=ForCowh$trt, mean=ForCowh.mean$fit,
   se=ForCowh.mean$se.fit )
ForCowh.mean[ !duplicated( ForCowh$trt ), ]
# least squares means
lsmean( ForCowh.aov )

# polynomial contrasts
coef( ForCowh.aov ) / sqrt(c(1,10,14,10,70))	# or coef( ForCowh.lm )
# use summary.lm to get standard errors
summary.lm( ForCowh.aov )		# or summary( ForCowh.lm )

# S-Plus plot of anova object
# plot( ForCowh.aov )

# residual plot with symbols
tmpdata = data.frame( x = predict( ForCowh.aov ),
  y = resid( ForCowh.aov ), group=ForCowh$trt )
print( xyplot( y ~ x, tmpdata, groups = group,
  xlab="mean by trt", ylab="residual" ))
# residual plot with jittered symbols
# add horizontal lines for 0 and +/- 1 SD
tmpdata$x = jitter( tmpdata$x )
print( xyplot( y ~ x, tmpdata, groups = group,
  panel = function(...) {
    panel.superpose(...)
    panel.abline( h = 0, lty = 2 )
    panel.abline( h=c(-1,1) * std.dev( ForCowh.aov ), lty=3 )
  },
  xlab="mean by trt", ylab="residual" ))

# Full set of cow measurements (3 periods per cow )
data( Forage )

Forage$code <- factor( Forage$code )
Forage$trt <- ordered( Forage$trt )
Forage$per <- ordered( Forage$per, c("Early","Middle","Late") )
Forage$ldmi <- log10( Forage$dmi )

Forage.fit <- aov( ldmi~hc*trt*per+Error( cow ),
   Forage, qr = TRUE )
summary( Forage.fit )

# projections
Forage.proj <- proj( Forage.fit )
Forage.wpint <- apply( Forage.proj$cow[,
   c("hc","trt","hc:trt") ], 1, sum )
Forage.spint <- apply( Forage.proj$Within[,
   c("per","hc:per","trt:per","hc:trt:per") ], 1, sum )
Forage.mean <- Forage.proj[[1]][1]

# H:23.1 Forage Treatment by Period (Heifer or Cow)

# Separate Analysis by Period

ylabs <- c("dry matter intake (DMI)","","")
xlabs <- c("(a) early period","(b) middle period","(c) late period")
mains = c("","Figure H:23.1","")
yaxis <- seq( 12, 24, 2 )
ylim <- log10( c(11.5,24.5) )
lper = levels( Forage$per )

for ( i in seq( along = lper )) {
   tmpfor <- Forage[Forage$per== lper[i], ]
   fit <- lm( log10( dmi )~trt*hc, tmpfor )
   ci.plot( fit, type="b",
     panel = function(x,y,...) {
       panel.abline( h=Forage.mean, lty=3 )
       panel.superpose(x,y,...)
     },
     xlab=xlabs[i], ylab=ylabs[i], main = mains[i],
     more = ( i != 3 ), split=c(i,1,3,1))
} 
#   axis( 2, log10( yaxis ), yaxis )
} 
% Converted by  version .
\keyword{datasets}
