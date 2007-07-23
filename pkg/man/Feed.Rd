\name{Feed}
\alias{Feed}
\title{Feed Intake and Weight Gain}
\usage{
data(Feed)
}
\format{
    Feed data frame with 12 observations on 7 variables.
    \tabular{rlll}{
	[,1] \tab trt \tab factor  \tab treatment identifier \cr
	[,2] \tab bwt \tab numeric \tab body weight \cr
	[,3] \tab bwg \tab numeric \tab body weight gain \cr
	[,4] \tab fi  \tab numeric \tab feed intake \cr
	[,5] \tab fe  \tab numeric \tab feed(?) \cr
	[,6] \tab cla \tab numeric \tab CLA amount \cr
 	[,6] \tab la  \tab numeric \tab LA amount
    }
}
\description{A scientist in food microbiology and toxicology (Chin et
al. 1994) examined the effect of additives on weight gain in animals.
He considered two levels (0.25\% and 0.5\%) of an additive \code{CLA},
a control and another additive (0.5\% \code{LA}), which were labelled,
respectively, as 1, 2, 3, 4.  The initial concern was the tremendous
differences in spread of \code{weight gain} for the four treatments.
Plots of \code{weight gain} against \code{feed intake} (Figure
\ref{fig:feedcov}) showed that this could be readily explained by
\code{feed intake}.  The \code{feed intake} slope is significant
(\emph{p} = 0.0088), as is the difference between \code{LA} and
\code{CLA} (\emph{p} = 0.0082) and the linear effect of \code{\%CLA}
(\emph{p} = 0.013). The \emph{p}-values are all computed using a Type
III approach.  Further, the effect of additive level seems to be linear
and uncorrelated with the feed intake.}
\source{SF Chin}
\references{Chin SF, Storkson JM, Albright KJ, Cook ME and Pariza MW (1994)
	`Conjugated linoleic acid is a growth factor for rats as
	shown by enhanced weight gain and imporved feed efficiency',
	\emph{J. Nutrition 124}, 2344-2349.}
\examples{
data( Feed )
Feed$trt <- factor( Feed$trt )
levels( Feed$trt ) = c("25CLA","5CLA","CTRT","LA5")

# adjusted treatment -- analysis of covariance
Feed.ancova <- aov( bwg ~ trt + fi, Feed )
# regression type fit
Feed.reg <- aov( bwg ~ fi + la + cla, Feed )
# simple regression
Feed.fi <- aov( bwg ~ fi, Feed)

# removing factor or covariate
Feed$bwgfi <- resid( Feed.fi ) + lsmean( Feed.fi )$pred
Feed.fi2 <- aov( bwgfi ~ fi + trt, Feed )
Feed.fi1 <- aov( bwgfi ~ trt, Feed )
Feed.trt <- aov( bwg ~ trt, Feed )
Feed$bwgtrt <- resid( Feed.trt ) + lsmean( Feed.trt, fac = NULL )$pred
Feed.trt2 <- aov( bwgtrt ~ fi + trt, Feed )
Feed.trt1 <- aov( bwgtrt ~ fi, Feed )

# checking for any interaction
Feed.int <- aov( bwg ~ fi * trt, Feed )
Feed$res <- resid( Feed.ancova )
Feed.res <- aov( res ~ fi * trt, Feed )
# checking for high interaction
Feed$hi <- as.numeric((Feed$cla + Feed$la) == .5)
Feed.hi <- aov( bwg ~ fi + trt + fi:hi, Feed )

# Figure F:17.1 Feed treatment adjusted for intake
## regression lines in following plot are not parallel
## want reg lines from Feed.reg
print( xyplot( bwg ~ fi, Feed, groups = trt,
  type = "p", pch = levels( Feed$trt ), cex = 2,
  col = 1 + seq( length( levels( Feed$trt ))),
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    pred = predict( Feed.reg )
    pch = levels( Feed$trt )
    col = 1 + seq( along = pch )
    lty = c(3,4,1,2)
    for( i in col - 1 ) {
      ii = Feed$trt == pch[i]
      panel.lines( x[ii], pred[ii], col = col[i], lty = lty[i] )
    }
    tmp = mean( x )
    panel.abline( v = tmp, lty = 4 )
    # least squares means at same
    panel.points( rep( tmp, 4 ), predict( Feed.reg,
      data.frame( la = Feed$la[!duplicated(Feed$trt)],
                 cla = Feed$cla[!duplicated(Feed$trt)],
                 fi = rep( tmp, 4 ) ) ), pch = 18 )
    panel.points( tapply( x, Feed$trt, mean ),
      tapply( y, Feed$trt, mean ), pch = 0 )
  },
  xlab = "Feed intake", ylab = "weight gain",
  main = "Figure F:17.1 Feed treatment adjusted for intake" ))
data.frame(se.bar( 2700, 1350, std.dev( Feed.ancova ), cap = "SD" ))

# Figure F:17.3 Feed ancove residuals with plot symbols
tmpdata = data.frame( x = fitted( Feed.ancova ),
  y = resid( Feed.ancova ), trt = Feed$trt )
print( xyplot( y ~ x, tmpdata, groups = trt, type = "p",
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    panel.abline( h = 0, lty = 2 )
    panel.abline( h = c(-1,1) * std.dev( Feed.ancova ), lty = 3 )
  },
  xlab = "fitted", ylab = "residual",
  main = "Figure F:17.1 Feed treatment adjusted for intake" ))

tmpfn = function(x,y,lty=c(3,4,1,2), col = 2:5 ) {
  pch = levels( Feed$trt )
  for( i in seq( along = pch )) {
    ii = pch[i] == Feed$trt
    panel.lines( x[ii], y[ii], lty = lty[i], col = col[i] )
  }
}
# Figure F:17.4 Feed residuals and confounding
print( xyplot( bwgfi ~ fi, Feed, group = trt, type = "p",
  col = 2:5, pch = levels( Feed$trt ), cex = 2,
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    tmpfn( x, fitted( Feed.fi1 ), rep(1,4) )
    tmpfn( x, fitted( Feed.fi2 ), rep(2,4) )
    tmp <- mean( x )
    ## vertical line at fi mean
    panel.abline( h = lsmean(Feed.fi)$pred, lty = 3 )
    panel.abline( v = tmp, lty = 3 )
    ## least squares means at same
    panel.points( rep( tmp, 4 ), predict( Feed.fi2,
        data.frame( trt = unique( Feed$trt ), fi = rep( tmp, 4 ))),
      pch = 18 )
    panel.points( tapply( x, Feed$trt, mean ),
      tapply( y, Feed$trt, mean ), pch = 0 )
  },
  xlab = "(a) weight gain adjusted for Feed intake",
  ylab = "residual weight gain",
  main = "Figure F:17.4" ),
  more = TRUE, split = c(1,1,2,1) )
data.frame(se.bar( 2700, 1475, std.dev( Feed.fi2 ), cap = "SD" ))

print( xyplot( bwgtrt ~ fi, Feed, group = trt, type = "p",
  col = 2:5, pch = levels( Feed$trt ), cex = 2,
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    tmpfn( x, fitted( Feed.trt1 ), rep(1,4), rep(1,4) )
    tmpfn( x, fitted( Feed.trt2 ), rep(2,4) )
    tmp <- mean( x )
    panel.abline( h = lsmean(Feed.trt)$pred, lty = 3 )
    panel.abline( v = tmp, lty = 3 )
    panel.points( rep( tmp, 4 ), predict( Feed.trt2,
        data.frame( trt = unique( Feed$trt ), fi = rep( tmp, 4 ))),
      pch = 18 )
    panel.points( tapply( x, Feed$trt, mean ),
      tapply( y, Feed$trt, mean ), pch = 0 )
  },
  xlab = "(b) weight gain adjusted for treatment",
  ylab = "residual weight gain",
  main = "Feed residuals and confounding" ),
  split = c(2,1,2,1) )
data.frame(se.bar( 2700, 1350, std.dev( Feed.trt2 ), cap = "SD" ))

# Figure F:17.5 Feed evidence for different slopes
print( xyplot( bwg ~ fi, Feed, group = trt, type = "p",
  col = 2:5, pch = levels( Feed$trt ), cex = 2,
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    tmp <- mean( x )
    panel.abline( v = tmp, lty = 3 )
    panel.points( rep( tmp, 4 ), predict( Feed.int,
        data.frame( trt = unique( Feed$trt ), fi = rep( tmp, 4 ))),
      pch = 18 )
    panel.points( tapply( x, Feed$trt, mean ),
      tapply( y, Feed$trt, mean ), pch = 0 )
    tmpfn( x, fitted( Feed.int ))
  },
   xlab = "(a) different slopes", ylab = "weight gain",
  main = "Figure F:17.5" ),
  more = TRUE, split = c(1,1,2,1) )
data.frame(se.bar( 2700, 1325, std.dev( Feed.int ), cap = "SD int" ))
data.frame(se.bar( 2690, 1325, std.dev( Feed.ancova ), cap = "add SD", adj = 1 ))

Feed$pred = fitted( Feed.ancova )
print( xyplot( res ~ pred, Feed, group = trt, type = "p",
  col = 2:5, pch = levels( Feed$trt ), cex = 2,
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    tmp <- mean( x )
    panel.abline( h = 0, lty = 5 )
    panel.points( tapply( x, Feed$trt, mean ),
      tapply( y, Feed$trt, mean ), pch = 0 )
    tmpfn( x, fitted( Feed.res ))
  },
  xlab = "(b) pure interaction", ylab = "additive residual",
  main = "Feed evidence for different slopes" ),
  split = c(2,1,2,1) )
} 
% Converted by  version .
\keyword{datasets}
