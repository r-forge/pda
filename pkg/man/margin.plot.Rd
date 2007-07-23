\name{margin.plot}
\alias{margin.plot}
\alias{margin.plot.default}
\title{Margin Plots of Response for Pairs of Factors}
\usage{
margin.plot(object, \dots)
\method{margin.plot}{default}(object, response, group, orient = "default", average = "full",
    xlab = xlabs, ylab = ylabs, \dots)
}
\description{}
\arguments{
  \item{object}{factor whose margins are to be plotted on the
    horizontal axis. Factor levels are plotted below the marginal means
    or effects.}
  \item{response}{numeric vector.}
  \item{group}{factor whose levels
    will be separate traces.}
  \item{\dots}{optional specification of other parameters to be passed
    to plotting routines.}
  \item{orient}{\code{"default"} = default (as in interaction.plot),
    \code{"switch"} = switch plot orientation of \code{object} and
    \code{group},
    \code{"both"} = plot both ways,
    \code{"full"} = use response,
    \code{"reduced"} = use predicted values
    from reduced (additive) model.
  }
  \item{xlab,ylab}{x and y axis labels}
  \item{average}{model to use for marginal averages}
}
\value{a plot will be created showing the requested function of
  responses for each level of the \code{object} at each level of the
    \code{group}. By default, lines for each value of the \code{group}
    are drawn in different styles so that they may be more easily
    distinguished.}
\seealso{\code{\link{tukey.plot}}, \code{\link{margin.plot.lm}},
    \code{\link{Tukey}}.}
\examples{
data( Tukey )
attach( Tukey )
margin.plot( lat, y, geno, orient = "switch", error = "red",
   xlab = "(b) geno marginals by lat", ylab = "",
   xpos = 115, ypos = 100, mod = "")
detach()}
\keyword{design}
% Converted by  version .
