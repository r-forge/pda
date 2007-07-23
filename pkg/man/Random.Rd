\name{Random}
\alias{Random}
\title{Random Data Examples}
\usage{
data( Random )
}
\description{Consider a population of \emph{classes}, say the counties of Wisconsin,
with an arbitrary class in this population labeled by \emph{u}.  This class
represents a sub-population of \emph{elements}, say the population of
farms in county \emph{u}.  If this class is selected for the experiment, an
independently selected random sample of its elements (farms) would
be selected for measurement, for instance of annual milk production per
farm.  The key questions concern estimating the average response
across the population (state) and understanding the sources of
variation in this measurement.  Note that a different way to run the
experiment would be to take a random sample of all elements (all farms
in the state), without regard to which class (county) they come from.
This latter design structure is completely
randomized, while the
former has two stages of randomization.  Substantial variation among
classes (counties) would be confounded with element-to-element
variation (among farms) in the latter design.
}
\seealso{\code{\link{Rantwo}}}
\examples{
data( Random )

# Figure G:19.1 Random model histograms
tmpdata = list( "(c) centered element values" = mean( Random ) +
    unlist( apply( Random, 2, function(x){x-mean(x)} )),
  "(b) class means" = apply( Random, 2, mean ),
  "(a) responses" = unlist( Random ))
tmpg = NULL
for( i in names( tmpdata ))
  tmpg = c( tmpg, rep( i, length( tmpdata[[i]] )))
tmpdata = data.frame( x = unlist( tmpdata ), g = tmpg )
histogram( ~ x | g, tmpdata, layout=c(1,3), nint = 30 )

# Figure G:19.2 Random model boxplots and tree diagram
tmp <- order( apply( Random, 2, median ) )
boxplot( as.list( Random[ , tmp ] ) )
#axis( 1, tmp, names( Random ) )
mtext( "(a) box-plots by class", 1, 2 )

plot( c(0.5,3.5), range( Random ), type = "n",
   xaxt = "n", xlab = "", ylab = "" )
mtext( c("mean","class","element"), 1, 0.5, at = 1:3 )
mtext( "(b) tree diagram", 1, 2 )
mtext( "response", 2, 2 )

Random.mean <- mean(unlist(Random))
Random.class <- apply(Random,2,mean)
points( 1, Random.mean )
points( rep( 2, length( Random ) ), Random.class )
points( rep( 3, length( unlist( Random ) ) ), unlist( Random ) )
branch( Random.mean, Random.class, 1:2, flip = TRUE )
for ( i in names( Random ) )
   branch( mean( Random[[i]] ), Random[[i]], 2:3, flip = TRUE )

# set up data for formal analysis
Random.frame <- data.frame( response = c( as.matrix( Random )))
Random.frame$class <- factor( rep( 1:10, rep( 10,10 )))

# one-factor anova (using aov() and fixed effects )
summary( aov( response ~ class, Random.frame ))
# one-factor anova (using aov() and random effects )
Random.ranaov <- aov( response ~ Error(class), Random.frame )

# one-factor anova (using lme() and random effects )
library( lme4 )
Random.lme <- lmer(response ~ 1 + (1|class), Random.frame )
summary( Random.lme )
VarCorr( Random.lme )
} 
% Converted by  version .
\keyword{datasets}

