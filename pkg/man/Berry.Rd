\name{Berry}
\alias{Berry}
\title{Cranberry Energy Balance}
\usage{
data(Berry)
}
\format{
  Berry data frame with 84 observations on 4 variables.
  \tabular{rlll}{
    [,1] \tab sample \tab factor \tab field sample or block \cr
    [,2] \tab tissue \tab factor \tab tissue type (flowering or vegetative) \cr
    [,3] \tab conc \tab numeric \tab starch concentration\cr
    [,4] \tab date \tab factor \tab Julian date for sample\cr
  }
}
\description{
This problem comes from a horticulturist
(Hagidimitriou and Roper 1994)
investigating the relative energy balance between \code{veg}etative and
\code{flo}wering parts of cranberry uprights in the marshes of northern
Wisconsin.  It is very difficult to separate individual plants as they
tend to overlap considerably.  In fact, a large cranberry bed may be
essentially one plant!

The researcher examined one cranberry bed which was divided into six
\code{sample}s for measurement purposes.  Soluble sugar and starch
(\code{conc}) were measured at various \code{period}s through the
growing season for random samples of roots, stems and uprights.  Later
in the season, the upright portions were distinguished as either
flowering or vegetative.  For each measurement day, the scientist
wished to compare photosynthesis in these plant parts.  The scientist
was interested in assessing the change in sugars and starches over the
course of the summer.  That is, was one consistently higher than the
other, or does this depend on \code{Julian} date?
}
\source{M Hagidimitriou}
\references{
  Hagidimitriou M and Roper T (1994)
  ``Seasonal changes in nonstructural carbohydrates in
  cranberry'',
  \emph{J. American Society of Horticultural Science 119},
  1029-1033.
}
\examples{
data( Berry )
library( lme4 )
Berry.lme <- lmer( conc ~ tissue * date + (1|sample),
   data = Berry, method = "REML" )
print( Berry.lme )
anova( Berry.lme )
}
% Converted by  version .
\keyword{datasets}
