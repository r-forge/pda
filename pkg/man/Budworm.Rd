\name{Budworm}
\alias{Budworm}
\title{Counting Data: Sex Response to Dose}
\usage{
data(Budworm)
}
\format{
    A data frame with 194 observations on 4 variables.
    \tabular{rlll}{
	[,1] \tab ldose \tab   numeric \tab log base 2 of dose\cr
	[,2] \tab dead \tab    numeric \tab number dead of 20\cr
	[,3] \tab sex \tab     factor \tab sex (M/F)
}}
\description{An experiment examined the toxicity of the pyrethroid
    \emph{trans}=cypermethrin on the tobacco budworm \emph{Heliothis
	virescens}. Batches of 20 budworm moths of each sex were exposed
    to the toxin for three days; the number dead or knocked down were
    reported.}
\source{Venables and Ripley (1999) \emph{Modern Applied Statistics
	with S-PLUS}, 3rd ed., ch. 7.}
\references{Collett (1991, p. 75) \emph{Modelling Binary Data},
    Chapman \& Hall, London.}
\examples{
data( Budworm )

Budworm.fit <- glm( dead ~ sex * ldose, data = Budworm,
   weight = total, family = binomial )
summary( Budworm.fit )

# Plot of Logit by Sex

plot( c(1,32), c(0,1), type = "n", xlab = "dose",
   ylab = "prob", log = "x" )
text( 2^Budworm$ldose, Budworm$dead, as.character( Budworm$sex ))
Budworm.ld <- seq( 0, 5, 0.1 )
for( i in levels( Budworm$sex ))
   lines( 2^Budworm.ld, predict( Budworm.fit,
      data.frame( ldose = Budworm.ld,
         sex = factor( rep (i, length( Budworm.ld )),
         levels = levels( Budworm$sex ))),
      type = "response" ))
rm( Budworm.ld )

# Parameters for Lines

summary( glm( dead ~ sex + ldose - 1, data = Budworm,
   weight = total, family = binomial ))$coefficients
} 
% Converted by  version .
\keyword{datasets}
