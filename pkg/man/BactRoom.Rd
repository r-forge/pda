\name{BactRoom}
\alias{BactRoom}
\title{Bacteria Virus in Plant Sap by Room}
\usage{
data(BactRoom)
}
\format{
  BactRoom data frame with 8 observations on 7 variables.
  \tabular{rlll}{
    [,1] \tab bact \tab  factor \tab  bacteria or control \cr
    [,2] \tab temp \tab  factor \tab  room temperature \cr
    [,3] \tab inoc \tab  factor \tab when innoculated\cr
    [,4] \tab bill \tab  numeric \tab bill length\cr
    [,5] \tab leg \tab   numeric \tab leg lenght\cr
    [,6] \tab code \tab  factor \tab plot code\cr
    [,7] \tab trt \tab   trt \tab bact*trt combination
  }
}
\description{
  An experiment performed at the National Wildlife Health Research
  Center under the supervision of Michael Samuel (Samuel, Goldberg, Thomas and
  Sharp 1995) examined the effect of a certain bacteria
  strain (\code{mycoplasma}) on the development of birds.  Unfortunately,
  due to the danger of aerial infection, they had to isolate treatment
  groups in different rooms.  A very conservative approach would take
  each chamber as an experimental unit, regardless of the number of
  birds per room.  There were four treatments: \code{cold, cold+myco,
    warm} and \code{warm+myco}.  In addition, the experiment was run in two
  \code{runs} since there were only \code{4} rooms.

  Measurements were taken on
  almost \code{200} chicks in these rooms.  It is possible to think of
  the experiment as having chicks (or eggs) randomly assigned
  to the rooms, and to consider the sample in each room as a random
  sample from a population of chicks exposed to that environment
  (combination of \code{temperature} and presence or absence of
  \code{mycoplasma} bacteria).
  This experiment on bird development was conducted in two \code{runs}
  separated by several weeks.  Several things could have changed in that
  time, including the \code{mycoplasma} culture, seasonal changes of
  chick growth and food or water conditions.  The scientist
  \code{inoculated} eggs in the first \code{run}, but decided to switch
  to inoculating young \code{chicks} in the later \code{run}.  Earlier
  analysis in the text has assumed that \code{inoculation} could be just
  considered as another factor.  Here it is viewed as a blocking factor
  with no replication.  That is, strictly speaking it is not possible to
  assess the main effect of \code{inoculation} method since there is no
  replication of \code{runs}.  However, it would be possible to assume
  that interactions with \code{run} were interactions with
  \code{inoculation}.
}
\source{Michael D Samuel (\url{mailto:michael_samuel@nbs.gov}),
  National Wildlife Health Center, Madison, WI
  (\url{http://www.emtc.nbs.gov/nwhchome.html})
}
\references{
  Samuel MD, Goldberg DR, Thomas CB and Sharp P (1995)
  ``Effects of \emph{Mycoplasma anatis} and cold stress on
  hatching success and growth of mallard ducklings,''
  \emph{J of Wildlife Diseases 31}, 172-178.
}
\seealso{\code{\link{Bacteria}}}
\examples{
data( BactRoom )
BactRoom$trt <- factor( BactRoom$trt )
BactRoom$code <- factor( BactRoom$code )
contrasts( BactRoom$temp ) <- contr.sum
contrasts( BactRoom$bact ) <- contr.sum
contrasts( BactRoom$inoc ) <- contr.sum

BactRoom.badd <- aov( bill ~ temp + bact + inoc, BactRoom )
summary( BactRoom.badd )

BactRoom.full <- aov( bill ~ temp * bact, BactRoom )
lsmean( BactRoom.full )
summary( BactRoom.full )

BactRoom.bill <- aov( bill ~ temp * bact * inoc, BactRoom )
lsmean( BactRoom.bill )

# C:8.1 Bacteria Interaction Plots

BactRoom.trt <- aov( bill ~ trt * inoc, BactRoom )
int.plot( BactRoom.trt, split = c(1,1,2,1), more = TRUE,
  bar.plot="none",
  xlab="(a) treatment by inoc", ylab="bill length (mm)",
  main = "Figure C:8.1a Wrong Way" )

int.plot( BactRoom.bill, bar.plot="none", split = c(2,1,2,1),
  xlab="(b) temp by bact", ylab="", lty = c(1,1,2,2),
  main = "Figure C:8.1b Right Way", white = FALSE )

# C:8.2 Bacteria Two Factor Interaction Plots

lsd.plot( BactRoom.full,
  xlab="(a) interaction plot", ylab="bill length (mm)",
  ypos=40, lty = 1,
  main = "Figure C:8.2a Interaction with LSD",
  split = c(1,1,2,1), more = TRUE )

margin.plot( BactRoom.full, ylim=range( BactRoom$bill ),
  effects=TRUE, xlab="(b) margin plot", ylab="",
  ypos=40, lty = 1, split = c(2,1,2,1),
  main = "Figure C:8.2b Margin Plot with LSD" )

# C:8.3 Bacteria Three Factor Interaction Plots

Bacteria.ylab <- c("","bill length (mm)")
Bacteria.xlab <- c("(a)","(b)")
Bacteria.main = c( "Figure C:8.3", "Three-Factor Interactions" )
Bacteria.inoc = levels( BactRoom$inoc )
for ( i in 1:2) {
  tmpdata = BactRoom[ BactRoom$inoc == Bacteria.inoc[i], ]
  print( xyplot( bill ~ temp, tmpdata, groups = bact, type="b",
    xlab = paste( Bacteria.xlab[i], Bacteria.inoc[i], "inoc" ),
    ylab = Bacteria.ylab[i],
    main = Bacteria.main[i] ),
    split = c(i,1,2,1), more = ( i == 1 ))
}
rm( Bacteria.xlab, Bacteria.ylab, Bacteria.main, Bacteria.inoc, tmpdata )

# C: Bacteria Half-Normal and Effect Plots
# order -- temp:bact:inoc temp:bact bact:inoc temp:inoc inoc bact temp

BactRoom.coef = sort( abs( coef( BactRoom.bill )[-1] ))
BactRoom.coef = BactRoom.coef * sqrt( length( BactRoom.coef ))

qhalfnorm = function(p)qnorm((1+p)/2)
qhalfquant <- qhalfnorm(ppoints(length(BactRoom.coef)))
print( qqmath( ~ BactRoom.coef, distribution = qhalfnorm,
  prepanel = prepanel.qqmathline,
  xlim = range(qhalfquant) + c(-.2,.2), ## make room for text
  panel = function(x, ...) {
    panel.qqmathline(x, distribution = qhalfnorm, lty = 2 )
    panel.text(qhalfquant, x, names(BactRoom.coef))
  },
  xlab = "(a) half-normal quantiles", ylab = "effects",
  main = "Figure C:8.4a Half-Normal Plot" ))

effect.plot( BactRoom.bill, xlim = c(.5,7.5), xaxt = "n",
   xlab = "(b) terms", ylab = "MS adjusted effects" )
title( "Figure C:8.4b Effects Plot" )
axis( 1, seq(1,7,by=2), c("T","I","IT","BIT") )
axis( 1, seq(2,6,by=2), c("B","BT","BI") )
}
\keyword{datasets}
