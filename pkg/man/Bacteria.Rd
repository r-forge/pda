\name{Bacteria}
\alias{Bacteria}
\title{Bacteria virus in plant sap}
\usage{
data(Bacteria)
}
\format{
    Bacteria data frame with 192 observations on 8 variables.
    \tabular{rlll}{
	[,1] \tab id \tab  factor \tab  bird identifier \cr
	[,2] \tab bact \tab  factor \tab  bacteria or control \cr
	[,3] \tab temp \tab  factor \tab  room temperature \cr
	[,4] \tab inoc \tab  factor \tab when innoculated\cr
	[,5] \tab bill \tab  numeric \tab bill length\cr
	[,6] \tab leg \tab   numeric \tab leg length\cr
	[,7] \tab code \tab  factor \tab plot code\cr
	[,8] \tab trt \tab   trt \tab bact*trt combination
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
<P>
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
\code{inoculation}.}
\source{Michael D Samuel (\url{mailto:michael_samuel@nbs.gov}),
    National Wildlife Health Center, Madison, WI
    (\url{http://www.emtc.nbs.gov/nwhchome.html})}
\references{Samuel MD, Goldberg DR, Thomas CB and Sharp P (1995)
	``Effects of \emph{Mycoplasma anatis} and cold stress on
	hatching success and growth of mallard ducklings,''
	\emph{J of Wildlife Diseases 31}, 172-178.}
    \seealso{\code{\link{Bacteria}}}
\examples{
# full data set--code not ready yet!
# model fits for bill and leg
data( Bacteria )
Bacteria$id <- factor( Bacteria$id )
Bacteria$cage <- factor( Bacteria$cage )
Bacteria$day <- ordered( Bacteria$day )
contrasts( Bacteria$temp ) <- contr.sum
contrasts( Bacteria$bact ) <- contr.sum
contrasts( Bacteria$inoc ) <- contr.sum

Bacteria.bill <- aov( bill ~ temp * bact * inoc, Bacteria )
Bacteria.leg <- aov( leg ~ temp * bact * inoc, Bacteria )

Bacteria.rcbd <- aov( bill ~ inoc + bact * temp, Bacteria )

# D:10.2 Bacteria Effect Plot by Bird

effect.plot( Bacteria.bill, xaxt = "n" )
title( "D:10.2 Bacteria Effect Plot by Bird" )
axis( 1, 1:8, c("temp","bact","inoc","T:B","T:I","B:I","T:B:I","resid") )

# H:22.2 Bacteria Interaction Plot Adjusted for Runs

effect.plot( Bacteria.rcbd,
   xlab="(a) terms", ylab="MS adjusted effects",
   xlim = c(0.5,5.5), xaxt = "n",
   width = .2 )
title( "H:22.2a Effect Plot Adjusted by Run" )
axis( 1, 1:5, c("run","bact","temp","inter","resid") )

margin.plot( Bacteria.rcbd, factors = c("temp","bact"),
  xlab = "(b) temperature by bacteria",
  ylab = "bill length (mm)",
  main = "H:22.2a Margin Plot Adjusted by Run" )
}
% Converted by  version .
\keyword{datasets}

