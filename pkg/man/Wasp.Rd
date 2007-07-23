\name{Wasp}
\alias{Wasp}
\title{Wasp Shape of Queens and Workers}
\usage{
data(Wasp)
}
\format{
    Wasp data frame with 100 observations on 14 variables.
    \tabular{rlll}{
	[,1] \tab caste \tab  factor \tab Queen or Worker \cr
	[,2] \tab TL \tab numeric \tab thorax length \cr
	[,3] \tab WL \tab numeric \tab wing length \cr
	[,4] \tab HH \tab numeric \tab head height \cr
	[,5] \tab HW \tab numeric \tab head width \cr
	[,6] \tab TH \tab numeric \tab thorax height \cr
	[,7] \tab TW \tab numeric \tab thorax width \cr
	[,7] \tab G1L \tab numeric \tab gonadium lenght 1 \cr
	[,7] \tab G2Wa \tab numeric \tab gonadium width 2a \cr
	[,7] \tab G2L \tab numeric \tab gonadium length 2 \cr
	[,4] \tab HL \tab numeric \tab head length \cr
	[,7] \tab G1Wb \tab numeric \tab gonadium width 1b \cr
	[,7] \tab G1Wa \tab numeric \tab gonadium width 1a\cr
	[,7] \tab G1H \tab numeric \tab gonadium height 1 \cr
    }
}
\description{Scientists studying social insects such as wasps and ants have
noted that queen and worker \code{caste}s can have very different sizes
and shapes.  For some species, the queens are simply larger, suggesting
that they continue to grow with the same basic `plan' as workers.
They are just fed for longer.  Jeanne, Graf and Yandell (1995)
examined a wasp species which does not follow this pattern.  They
measured 13 responses on 50 workers and 50 queens.  The measurements
were coded based on body part -- head (\code{H}), thorax
(\code{T}), wing (\code{W}) or gonadium (\code{G}) -- and kind of
measurement -- width (\code{W}), height (\code{H}) and length (\code{L}).
The gonadium (also known as gastral tergite) has two length and three
width measurements.}
\source{Robert L Jeanne}
\references{Jeanne RL, Graf CA and Yandell BS (1995)
	`Non-size-based morphological castes in a social insect',
	\emph{Naturwissenschaften 82}, 296-298.}
\examples{
data( Wasp )

## Multivariate discriminant analysis
library(MASS)
Wasp.lda <- lda( Wasp[,-1], Wasp$caste )
Wasp.lda
Wasp$DA <- predict( Wasp.lda, Wasp[,-1] )$x

## Univariate t statistics
tmp <- Wasp$caste == "Q"
Wasp.t <- apply( Wasp[,-1], 2, function( x ) {
   ( mean( x[tmp] ) - mean( x[!tmp] ) ) /
   sqrt( ( var( x[tmp] ) + var( x[!tmp] ) ) / sum( tmp ) ) } )
rm( tmp )
Wasp.t <- sort( Wasp.t )
Wasp.tmpy <- cor( Wasp[,2:14], Wasp$DA )[ names( Wasp.t[-1] ), 1 ]

## Figure F:18.2 Wasp scatter plots

pairs( Wasp[ , c("HW","TW","G1Wa","DA") ], panel =
   function( x, y, ... ) {
     text( x, y, as.character( Wasp$caste ) ) },
   labels = c("head","thorax","gonadium",
     "discriminant\nanalysis" ) )
title( "Figure F:18.2 Wasp scatter plots" )

## Figure F:18.3 Wasp discriminant analysis

plot( c(1,-1) * Wasp.t["DA"], c(1,-1), type = "n",
   xlab ="t statistics", ylab = "correlation with DA" )
text( c(1,-1) * Wasp.t["DA"], c(1,-1), "DA" )
lines( Wasp.t[-1], Wasp.tmpy )
points( Wasp.t[-1], Wasp.tmpy )
text( Wasp.t[-1] + 5 * sign( Wasp.t[-1] ),
  uncollide(Wasp.tmpy,.05), names( Wasp.t[-1] ) )
abline(h=0,v=0,lty=2)
mtext("queen larger",1,2,adj=1)
mtext("worker larger",1,2,adj=0)
title( "Figure F:18.3 Wasp discriminant analysis" )

# Figure F:18.5 Wasp analysis of covariance

Wasp.cov <- aov( HW ~ G1Wa * caste, Wasp )
anova( Wasp.cov )
Wasp.par <- aov( HW ~ G1Wa + caste, Wasp )
anova( Wasp.par )

print(xyplot(HW ~ G1Wa, Wasp, group = caste, type = "p",
  pch = levels(Wasp$caste), col = 1 + seq(along = levels(Wasp$caste)),
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    lc = levels(Wasp$caste)
    for(i in seq( along = lc )) {
      ii = lc[i] == Wasp$caste
      panel.lines(Wasp$G1Wa[ii], predict(Wasp.cov)[ii],
        lty = 3, col = i)
      panel.lines(Wasp$G1Wa[ii], predict(Wasp.par)[ii],
        lty = 1, col = i)
    }
  },
  xlab ="(a) gonadium width (G1Wa)",
  ylab = "head width (HW)",
  main = "Figure F:18.5 Wasp" ),
  more = TRUE, split = c(1,1,2,1) )

## residual analysis after removing TW
Wasp.resid <- Wasp[,c("caste","HW","G1Wa","TW")]
Wasp.resid$Hres <- resid( aov( HW ~ TW, Wasp.resid ) )
Wasp.resid$Gres <- resid( aov( G1Wa ~ TW, Wasp.resid ) )
Wasp.covr <- aov( Hres ~ TW + G1Wa + caste, Wasp.resid )
anova( Wasp.covr )
Wasp.covwr <- aov( Hres ~ G1Wa + caste, Wasp.resid )
anova( Wasp.covwr )
Wasp.resid$TW <- rep( mean( Wasp.resid$TW ), nrow( Wasp.resid ) )

print(xyplot(Hres ~ G1Wa, Wasp.resid, group = caste, type = "p",
  pch = levels(Wasp.resid$caste), col = 1 + seq(along = levels(Wasp.resid$caste)),
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    lc = levels(Wasp.resid$caste)
    for( i in seq(along = lc)) {
      ii = lc[i] == Wasp.resid$caste
      panel.lines(Wasp.resid$G1Wa[ii], predict(Wasp.covr, Wasp.resid)[ii],
        lty = 1, col = i )
      panel.lines(Wasp$G1Wa[ii], predict(Wasp.covwr)[ii],
        lty = 3, col = i )
    }
  },
  xlab ="(b) gonadium width (G1Wa)",
  ylab = "HW residuals on TW",
  main = "Analysis of Covariance" ),
  split = c(2,1,2,1) )
} 
% Converted by  version .
\keyword{datasets}
