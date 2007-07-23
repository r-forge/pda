\name{Tomato}
\alias{Tomato}
\title{Tomato data to find genes for growth}
\usage{data(Tomato)}
\format{
    A data frame with 194 observations on 4 variables.
    \tabular{rlll}{
	[,1] \tab entry \tab   factor \tab  identifier for line entry\cr
	[,2] \tab yr \tab      factor \tab  year of measurement\cr
	[,3] \tab mfwlog \tab  numeric \tab log( mean flowering time )\cr
	[,4] \tab tg430 \tab   ordered \tab A<H<B genotype at marker TG430
}}
\description{A plant scientist is interested in finding out the location of
genes which control various plant growth attributes.  He is using
molecular markers to do this.  He thinks he has located a major gene
for fruit weight (\code{fw}) close to marker \code{tg430}.

The original experiment covered two years (\code{yr}), with 93 unique
plant entries (\code{entry}).  There were several replications each
year.  Unfortunately, the data now available to the scientist consists
of the mean fruit weight (\code{mfw}), averaged across the replicates
for each year.  [The raw data are in notebooks halfway around the
world!]  The scientist believes a $\log_{10}$ transformation ({\tt
mfwlog}) is reasonable, and has presented that data in that way.

Nevertheless there are still 2 years of data for most entries.  The
marker \code{tg430} can be used to classify entries into one of three
categories, 1 = parent \code{A}, 3 = parent \code{B}, 2 = \code{H}ybrid
of \code{A} and \code{B} (and \code{.} = missing marker value).  The
scientist is particularly interested in the `additive effect' (parent
\code{A} -- parent \code{B}) and the `dominance effect' (\code{H}ybrid --
mean of parents).  Note that if the dominance is zero, then the hybrid
would be halfway between the two parents.

}
\source{Professor Irwin Goldman (\url{mailto:igoldman@facstaff.wisc.edu}),
	Horticulture Department, UW-Madison}
\references{Goldman IL, Paran I, Zamir D (1995)
	`Quantitative trait locus anlaysis of a recombinant inbred
	line population derived from a \emph{Lycopersicon esculentum}
	x \emph{Lycopersicon cheesmanii} cross',
	\emph{Theoretical & Applied Genetics 90}, 925-932.}
\examples{
data( Tomato )

# make sure entry, tg430 and yr are all factors
Tomato$entry <- factor( Tomato$entry )
Tomato$tg430 <- ordered( Tomato$tg430, c("A","H","B") )
Tomato$yr <- ordered( Tomato$yr )
# reduce to complete data
Tomato1 <- Tomato[ !( is.na( Tomato$tg430 )
   | is.na( Tomato$mfwlog ))
   & Tomato$yr == 1, ]

# Figure 4.1: Histograms
histogram( ~ mfwlog | tg430, Tomato1, nint = 30, layout = c(1,3),
  main = "Figure B:4.1. Tomato Histograms by Group" )

# Figure 4.2: Box-Plots
bwplot( mfwlog ~ tg430, Tomato1,
  xlab = "Tomato Allele Type", ylab= "Log Flower Time",
  main = "Figure B:4.2. Tomato Box-Plots by Group" )

# Figure 5.1: Confidence Intervals

# fit one-factor anova
Tomato.aov <- aov( mfwlog ~ tg430, Tomato1 )
Tomato.aov
summary( Tomato.aov )
# least squares means ( uses library( pda ) )
lsmean( Tomato.aov )

# 95% confidence intervals by genotype ( uses library( pda ) )
ci.plot( Tomato.aov, level = 0.05,
   crit=qt( 1 - 0.05 / 2, df.resid( Tomato.aov )) / sqrt(2),
   xlab=paste("(a) ",100 * ( 1 - 0.05 ),
      "\% CI / sqrt(2)", sep = "" ),
   ylab = "log flower time",
   main = "Figure 5.1(a) Confidence Intervals",
  split = c(1,1,1,1) )

# notched box-plots to compare with CIs (NA in bwplot currently)
attach( Tomato1 )
boxplot( split( mfwlog, tg430 ), notch = TRUE,
   xlab = "(b) notched box-plots",
   ylab = "log flower time",
   main = "Figure 5.1(b) Tomato Notched Box-Plots" )
detach()
}
% Converted by  version .
\keyword{datasets}
