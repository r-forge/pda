\name{BrandX}
\alias{BrandX}
\title{Preference for Brand X Detergent}
\usage{
data(BrandX)
}
\description{These are results of a detergent brand preference study, with
    respondents classified by temperature and softness used during their
    wash cycle. Venables and Ripley call \code{Brand} the response
    factor and \code{M.user}, \code{Temp} and \code{Soft} stimulus
    factors. Thus they are most interested in \code{Brand} and in
    interactions that involve \code{Brand}.}
\source{Ries and Smith (1963) The use of chi-square for preference
    testing in multidimensional problems. \emph{Chemical Engineering
	Progress 59}, 39-43.}
    \references{
      \item{Venables and Ripley (1999)}{\emph{Modern Applied Statistics
	with S-PLUS}, 2nd ed., ch. 7. (this is \emph{not} in the 3rd ed)}
    \item{Cox and Snell (1989)}{\emph{Analysis of Binary Data}, Chapman \& Hall.}}
\examples{
# Four-way Contingency Table

data( BrandX )

# fit poisson model

BrandX.fit <- glm( Fr ~ M.user*Temp*Soft+Brand,
   family = poisson, data = BrandX )
anova( BrandX.fit, test = "Chisq" )
drop1( BrandX.fit, formula( BrandX.fit), test = "Chisq" )

BrandX.step <- step(BrandX.fit,
   list( lower = formula( BrandX.fit ), upper = ~.^3 ),
   scale = 1, trace = FALSE)
BrandX.step$anova
anova( BrandX.step, test = "Chisq" )
BrandX.mod <- glm( terms( Fr ~ M.user*Temp*Soft +
   Brand*M.user*Temp, keep.order = TRUE ),
   family = poisson, data = BrandX )
summary(BrandX.mod, correlation = FALSE, test = "Chisq" )
} 
% Converted by  version .
\keyword{datasets}
