\name{tukey.plot}
\alias{tukey.plot}
\alias{mandel.plot}
\title{Tukey/Mandel Interaction Plots of Response for Pairs of Factors}
\usage{
tukey.plot(x.factor, response, group, \dots, orient = "default",
  xlab, ylab)
mandel.plot(x.factor, response, group, \dots, orient = "default",
  xlab, ylab)
}
\description{}
\arguments{
  \item{x.factor}{factor whose margins are to be plotted on the
    horizontal axis. Factor levels are plotted below the marginal means
    or effects.}
  \item{response}{numeric vector.}
  \item{group}{factor whose levels will be separate traces.}
  \item{\dots}{optional specification of other parameters to be passed
    to plotting routines.}
  \item{orient}{
    "default" = default (as in interaction.plot),
    "switch" = switch plot orientation of x.factor and trace.factor,
    "both" = plot both ways.}
  \item{xlab,ylab}{x and y axis labels}
}
\value{a plot will be created showing the requested function of
    responses for each level of the \code{x.factor} at each level of the
    \code{trace.factor}. By default, lines for each value of the
    \code{trace.factor} are drawn in different styles so that they may
    be more easily distinguished.}
\seealso{\code{\link{margin.plot}},
    \code{\link{Tukey}}.}
\examples{
data( Tukey )
attach( Tukey )
tukey.plot( geno, y, lat,
   xlab = "(b) geno marginals by lat", ylab = "" )
detach()
}
\keyword{design}
% Converted by  version .
