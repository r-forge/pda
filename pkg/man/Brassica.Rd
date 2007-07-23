\name{Brassica}
\alias{Brassica}
\title{Brassica Response to Herbicide}
\usage{
data(Brassica)
}
\format{
    Brassica data frame with 144 observations on 10 variables.
    \tabular{rlll}{
	[,1] \tab crop    \tab  factor \tab  crop identifier (1-7) \cr
	[,2] \tab pyrrat  \tab  ordered \tab pyridate level (0-2) \cr
	[,3] \tab block   \tab  factor \tab  block (1-3) \cr
	[,4] \tab inj7dat \tab  numeric \tab injury after 7 days \cr
	[,5] \tab in14dat \tab  numeric \tab injury after 14 days \cr
	[,6] \tab yldkga  \tab  numeric \tab yield in kilograms \cr
	[,7] \tab pconyld \tab  numeric \tab percent control yield \cr
	[,8] \tab pcnblk  \tab  numeric \tab percent control in block yield\cr
	[,9] \tab year    \tab  ordered \tab year - 1990\cr
	[,10] \tab trt     \tab  factor \tab treatment code
    }
}
\description{The scientist examined the effects of rate of
application of a certain herbicide (pyridate) on a range of closely
related \code{crop}s.  The experiment was conducted over two
\code{year}s, with three \code{block}s each \code{year}.
The two \code{year}s were conducted in roughly the same location,
but there was separate randomization in each \code{year}. Three
herbicide rates (\code{pyrrat} = 0, 1, 2) were randomly assigned
to separate rows within each \code{crop} plot.  Recorded values
consist of \code{yield} per row for each \code{crop*pyrrat}
combination within each \code{block} for each \code{year}.
Since \code{crop} 2 was not used in the second \code{year}, it
is dropped from further consideration.}
\source{Alan Miller, U WI Horticulture
	(\url{http://www.hort.wisc.edu})}
\references{Miller A (1993)
	`Pyridate tolerance in \emph{Brassicaceae} crops',
	PhD Dissertation, Department of Horticulture,
	University of Wisconsin-Madison.}
\examples{
# full data set--code not ready yet!
# model fits for bill and leg
data( Brassica )

# drop Crop 2
Brassica <- data.frame(Brassica[ Brassica$crop != 2, ])
Brassica$crop <- factor( as.character( Brassica$crop ))

# Analysis based on lmer()
library(lme4)
Brassica1 <- data.frame( Brassica[ Brassica$year == 1, ] )
Brassica.lme <- lmer( yldkga ~ crop * pyrrat + (1|crop:block),
   data = Brassica1 )
anova( Brassica.lme )
VarCorr( Brassica.lme )

# Analysis based on aov() using Error() specification

Brassica.aov <- aov ( yldkga ~ crop * pyrrat +
   Error( block/crop ), Brassica1 )
summary( Brassica.aov )

# Analysis of full data set: both years

Brassica.full <- aov ( yldkga ~ crop * pyrrat * year +
   Error( year:block + crop:year:block ), Brassica )

Brassica.red <- aov ( yldkga ~ crop * year + pyrrat + crop:pyrrat +
   Error( year:block + crop:year:block ), Brassica)

# Analysis based on lmer()
Brassica1 <- data.frame( Brassica[ Brassica$year == 1, ] )
Brassica.lme <- lmer( yldkga ~ crop * pyrrat + (1|block:crop),
   data = Brassica1 )
Brassica$yrblk = interaction(Brassica[,c("year","block")])
Brassica$cropyrblk = interaction(Brassica[,c("crop","year","block")])
Brassica.fulle = lmer(yldkga ~ crop * pyrrat * year +
  (1|yrblk:cropyrblk), data = Brassica )
anova(Brassica.fulle)
Brassica.rede = lmer(yldkga ~ crop * year + pyrrat + crop:pyrrat +
  (1|yrblk:cropyrblk), data = Brassica )
anova(Brassica.rede)
anova(Brassica.rede,Brassica.fulle)

lsd.plot(Brassica.rede,factor=c("pyrrat","crop"),
   xlab = "pyr rate by crop",
   ylab = "mean crop yield" )

\dontrun{
# The following code is specialized to get plots for the book.
# It is not recommended as standard practice, but presented to
# show one torturous way to get LSMEANS for a nested design.
Brassica.red <- aov (yldkga ~ crop * year + pyrrat + crop:pyrrat +
   Error(year:block + crop:year:block), Brassica, qr = TRUE)
Brassica.proj <- proj ( Brassica.red )
Brassica.ref <- data.frame(Reference=c(Brassica.proj[[1]]))
for (i in names(Brassica.proj)[-1])
   Brassica.ref[[i]] <- apply( Brassica.proj[[i]], 1, sum )

# whole plot
Brassica.wp <- nested( Brassica.proj, Brassica, c("block","year","crop"), "yldkga",
   c(1,3), Brassica.ref )

Brassica.wpaov <- aov( yldkga ~ year*block + crop + crop:year, Brassica.wp,
   weight = Brassica.wp$weight)

# split plot
Brassica.sp <- nested( Brassica.proj, Brassica, c("block","year","crop","pyrrat"),
   "yldkga", c(1,4), Brassica.ref )

Brassica.spaov <- aov( yldkga ~ year*block*crop + pyrrat +
   crop:pyrrat, Brassica.sp )

margin.plot( Brassica.wpaov, factors = c("crop","year"),
   xlab = "(a) crop by year", ylab = "yield in 10000 kg" )

#Figure H:24.1 Brassicaica split split plot anova

lsd.plot( Brassica.spaov, factors = c("pyrrat","crop"),
   xpos = .6, ypos = 3.25,
   xlab = "(b) pyr rate by crop",
   ylab = "deviation from mean crop yield" )
}
}
% Converted by  version .
\keyword{datasets}
