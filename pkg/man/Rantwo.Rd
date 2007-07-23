\name{Rantwo}
\alias{Rantwo}
\title{Random Two-Factor Data Example}
\usage{
data( Rantwo )
}
\description{
  Consider a population of seed packets and a population of fields in
  Dane  County.  The experiment consists of taking a random sample of
  seed packets and a random sample of fields, then placing several seeds
  from each seed packet in each field.  Individual plants grown from
  seed are measured for yield at the end of the growing season.
}
\references{}
\examples{
data( Rantwo )
Rantwo$seed = factor( Rantwo$seed )
Rantwo$farm = factor( Rantwo$farm )
Rantwo$seedfarm = with( Rantwo, seed:farm)
Rantwo.fit <- aov( yield ~ seed * farm, Rantwo )
summary( Rantwo.fit )
library( lme4 )
Rantwo.lme <- lmer(yield ~ 1 + (1|seed) + (1|farm) + (1|seedfarm),
   data = Rantwo )
summary( Rantwo.lme )
Rantwo.lmep <- lmer(yield ~ 1 + (1|seed) + (1|farm),
   data = Rantwo )
anova(Rantwo.lme,Rantwo.lmep)
}
% Converted by  version .
\keyword{datasets}

