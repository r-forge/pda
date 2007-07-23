\name{Drying}
\alias{Drying}
\title{Water Temperature Effect on Cow Dryinging}
\usage{
data(Drying)
}
\format{
    Drying data frame with 261 observations on 6 variables.
    \tabular{rlll}{
	[,1] \tab cow \tab  factor \tab  cow identifier \cr
	[,2] \tab period \tab  factor \tab  period (1-3) \cr
	[,3] \tab sample \tab  ordered \tab sample number per cow \cr
	[,4] \tab drying \tab  factor \tab  drying treatment (1-3) \cr
	[,5] \tab dmrem \tab  numeric \tab dry matter remaining \cr
	[,6] \tab cprem \tab numeric  \tab crude protein remaining
    }
}
\description{A dairy scientist wished to examine
the effect of forage sample drying technique on ruminal disappearance
of dry matter and crude protein.  He gathered 30 uniquely
different forage (silage) \code{sample}s, roughly ten (9--11) per
\code{period}.  Each \code{sample} was divided each into sub-samples for
\code{drying} by three different methods (\code{microwave}, \code{freeze}
dryer and convection \code{oven}).  There were three Holstein
\code{cow}s available with canulas (openings into the stomach for 
experimental observation), but each \code{cow}'s rumen could only hold
about \code{150} grams of forage.  The scientist wanted to make sure
that each \code{sample} was exposed to all nine combinations of the
three \code{cows} and three \code{drying} techniques.  Therefore,
each \code{drying} sub-sample was further divided into 5 gram
portions which were randomly assigned to the \code{cows}.
The percentages of dry matter (\code{dm}) and crude protein (\code{cp})
remaining after canulation were measured, although some portions were
inadvertently destroyed.}
\source{Pat Hoffman}
\references{Drendel TR, Hoffman PC, Bringe AN and Syverud TD (1994)
	`The effect of a premilking teat disinfectant on somatic cell
	count and clinical mastitis',
	Technical Report R3598, Research Division, CALS, University of
	Wisconsin-Madison.}
\examples{
data( Drying )
Drying$cow <- factor( Drying$cow )
Drying$period <- factor( Drying$period )
Drying$sample <- ordered( Drying$sample )
Drying$drying <- factor( Drying$drying )

tmpar = par( mfrow = c(3,3) )
Drying.xlim <- range( Drying$cprem, na.rm = TRUE )
Drying.ylim <- range( Drying$dmrem, na.rm = TRUE )
for( j in levels( Drying$drying ))
  for( i in levels( Drying$period )) {
    Drying.ij <- i == Drying$period & j == Drying$drying
    plot( Drying$cprem[Drying.ij], Drying$dmrem[Drying.ij],
      xlim = Drying.xlim, ylim = Drying.ylim, type = "n",
      xlab = "cprem",
      ylab = "dmrem" )
    title( paste( "drying", j, "period", i ))
    text( Drying$cprem[Drying.ij], Drying$dmrem[Drying.ij],
      interaction( Drying$sample[Drying.ij],
        unclass(Drying$cow[Drying.ij] )))
  }
mtext( "Drying Problem (sample.cow)", outer = TRUE )
par( tmpar )
rm( Drying.ij, Drying.xlim, Drying.ylim, tmpar )
} 
% Converted by  version .
\keyword{datasets}
