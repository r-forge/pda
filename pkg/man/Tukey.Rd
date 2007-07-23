\name{Tukey}
\alias{Tukey}
\alias{Mandel}
\title{Tukey \& Mandel Interactions}
\usage{
data(Tukey)
}
\format{
  Tukey data frame with 9 observations on 3 variables.
  \tabular{rlll}{
    [,1] \tab geno \tab  factor \tab genotype \cr
    [,2] \tab lat  \tab  factor \tab latitude \cr
    [,3] \tab y    \tab  numeric \tab response 
  }
}
\description{A plant geneticist wanted to compare the growth characteristics of two
  inbred lines \emph{A} and \emph{B} with a hybrid combination \emph{H}.  The plants
  were grown at field stations at three latitudes, 25, 35, 45, and a
  number of measurements were taken.  Unfortunately, the records were
  burned and the scientist only has the mean values from each station.
  The data suggest that genotype and latitude are not
  additive.  All three lines appear to perform equally well at latitude
  25.  Line \emph{A} grows shorter while line \emph{B} and the hybrid \emph{H}
  Grow taller as latitude increases.  The hybrid shows evidence of
  `hybrid vigor', doing much better than either parent line at higher
  latitudes.}
\references{Milliken GA and Johnson DE (1989)
  \emph{Analysis of Messy Data vol. 2: Non-replicated Experiments}.
  Chapman & Hall, London.
}
\seealso{\code{\link{tukey.plot}}, \code{\link{margin.plot}}.}
\examples{
###NOTE: These data differ somewhat from the book figures.
data( Tukey )
Tukey$geno <- ordered( Tukey$geno, c("A","H","B") )
Tukey$lat <- ordered( Tukey$lat )

# Figure C:9.3  Interaction Plots for Tukey Data

trellis.par.set(theme=col.whitebg())
print(xyplot(y ~ lat, Tukey, groups = geno, cex = 2,
               pch = levels( Tukey$geno ), type = "b",
   xlab = "(a) latitude by geno", ylab = "height",
   main = "Figure C:9.3" ),
   more = TRUE, split = c(1,1,2,1) )
print(xyplot(y ~ geno, Tukey, groups = lat, cex = 2,
               pch = levels( Tukey$lat ), type = "b",
   xlab = "(b) geno by latitude", ylab = "height",
   main = "Tukey Interaction Plot" ),
   split = c(2,1,2,1) )

# Figure C:9.4 Tukey raw margin plots

tmpar <- par( mfrow = c(1,2) )

attach( Tukey )
margin.plot( lat, y, geno, xlim = c(95,105), error = "red",
   xlab = "(a) lat marginals by geno", ylab = "height",
   xpos = 103.5, ypos = 100, mod = "",
   main = "Figure C:9.4",
   more = TRUE, split = c(1,1,2,1) )
margin.plot( lat, y, geno, orient = "switch", error = "red",
   xlab = "(b) geno marginals by lat", ylab = "height",
   xpos = 115, ypos = 100, mod = "",
   main = "Raw Margin Plot",
   split = c(2,1,2,1) )
detach()

# Figure C:9.5 Tukey fitted margin plots

attach( Tukey )
Tukey.fit <- tukey.plot( lat, y, geno,
   xlab = "(a) lat marginals by geno", ylab = "height",
   main = "Figure C:9.5",
   more = TRUE, split = c(1,1,2,1) )
tukey.plot( geno, y, lat,
   xlab = "(b) geno marginals by lat", ylab = "height",
   main = "Fitted Margin Plot",
   split = c(2,1,2,1) )
detach()

Tukey.fit$lsm
summary.aov( Tukey.fit$reduced )
summary.aov( Tukey.fit$full )
}
% Converted by  version .
\keyword{datasets}
