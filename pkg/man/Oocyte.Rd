\name{Oocyte}
\alias{Oocyte}
\title{Oocyte Development}
\usage{
data(Oocyte)
}
\format{
  Oocyte data frame with 192 observations on 8 variables.
  \tabular{rlll}{
    [,1] \tab fors \tab  factor \tab  forskolin level \cr
    [,2] \tab h89 \tab  factor \tab  H89 level \cr
    [,3] \tab date \tab  factor \tab  date \cr
    [,4] \tab dup \tab  factor \tab duplicate \cr
    [,5] \tab toocytes \tab  numeric \tab toocytes\cr
    [,6] \tab gv \tab   numeric \tab gv\cr
    [,7] \tab mi \tab  numeric \tab mi
  }
}
\description{Rose-Hellekant and Bavister (1996) has been studying the chemistry of oocyte (egg)
  development.  In one experiment, she examined the effect of two
  agents, \code{fors}kolin and \code{h89}.
}
\source{Terry Rose-Hellekant}
\references{Rose-Hellekant TA and Bavister BD (1996)
  `Roles of protein kinase A and C in spontaneous maturation
  and in forskolin or 3-isobutyl-1-methylxanthine maintained
  meiotic arrest of bovine oocytes',
  \emph{Molecular Reproduction & Development 44}, 241-249.
}
\examples{
data( Oocyte )

Oocyte$asper <- 100 * asin( sqrt( Oocyte$gv / Oocyte$toocytes ) )

Oocyte$h89 <- ordered( Oocyte$h89 )
Oocyte$fors <- ordered( Oocyte$fors )
Oocyte$date <- ordered( Oocyte$date )

sample.size( Oocyte$h89, Oocyte$fors )

cat( "Type I sums of squares:\n" )
Oocyte.fit <- aov( asper ~ h89 * fors, Oocyte )
summary( Oocyte.fit )
summary( aov( asper ~ fors * h89, Oocyte))
cat( "Type III sums of squares:\n" )
drop1( Oocyte.fit, formula( Oocyte.fit ), test = "F" )

cat( "Additive model:\n" )
Oocyte.afit <- aov( asper ~ h89 + fors, Oocyte )
cat( "Type II sums of squares:\n" )
drop1( Oocyte.afit, formula( Oocyte.afit ), test = "F" )

# D:10.1 Oocyte additive and interactive models

lsd.plot( Oocyte.afit, factors = c("fors","h89"),
  xlab = "(a) additive fit", ylab = "fors",
  more = TRUE, split = c(1,1,2,2) )
margin.plot( Oocyte.afit, factors = c("fors","h89"),
  xlab = "(b) additive margin", ylab = "",
  more = TRUE, split = c(2,1,2,2) )
lsd.plot( Oocyte.fit, factors = c("fors","h89"),
  xlab = "(c) interactive fit", ylab = "fors",
  more = TRUE, split = c(1,2,2,2) )
margin.plot( Oocyte.fit, factors = c("fors","h89"),
   xlab = "(d) interactive margin", ylab = "",
  split = c(2,2,2,2) )

# D:10.2 Oocyte Effects Plot
effect.plot( Oocyte.fit,
  main = "Figure D:10.2 Oocyte Effects Plot" )

cat( "Cautionary tale on Type III sums of squares\n",
   "Be sure to have sum to zero contrasts!\n" )
# first alter contrasts for H89
contrasts( Oocyte$h89 )
contrasts( Oocyte$h89 ) <- contr.treatment
contrasts( Oocyte$h89 )
Oocyte.trt <- aov( asper ~ h89 * fors, Oocyte )
drop1( Oocyte.trt, formula( Oocyte.trt ), test = "F" )

# now alter contrasts for forskolin
contrasts( Oocyte$fors )
contrasts( Oocyte$fors ) <- contr.treatment
contrasts( Oocyte$fors )
Oocyte.trt <- aov( asper ~ h89 * fors, Oocyte)
drop1( Oocyte.trt, formula( Oocyte.trt ), test = "F" )

# resetting contrasts for your sanity
contrasts( Oocyte$h89 ) <- contr.poly
contrasts( Oocyte$fors ) <- contr.poly
}
% Converted by  version .
\keyword{datasets}
