\name{margin.plot.lm}
\alias{margin.plot.lm}
\title{Margin Plots of Response for Pairs of Factors in Linear Model}
\usage{
\method{margin.plot}{lm}(object, reduced, data, factors, lsm, margin,
   error = "full", bar.plot = "lsd",
   xlab, ylab, \dots)}
\arguments{
  \item{object}{full model object from linear model fit}
  \item{reduced}{reduced model object, for computation of marginal
    means or effects(default=object)}
  \item{data}{data frame in which to interpret the variables named
    in the full and reduced objects (determined from reduced fit if
    omitted)}
  \item{factors}{Character string of length 2 with names of x.factor
    and trace.factoras found in fit and data (default is first 2
    factors in fit).}
  \item{lsm}{Least squares means for plotting (default is fit from
    full model).}\item{margin}{Least squares marginal means for
    horizontal axis(default is fit from reduced model rounded to
    proper precision).}
  \item{error}{model to use for error calculation ("full" or
    "reduced").}
  \item{bar.plot}{error bars:Fisher's LSD ("lsd"=default)or confidence
    intervals ("ci")or none ("none").}
  \item{xlab,ylab}{plot parameters}
  \item{\dots}{optional specification of other parameters to be passed
    to plotting routines.}}
\value{a plot will be created showing the requested function of
  responses for each level of the \code{x.factor} at each level of the
  \code{trace.factor}. By default, lines for each value of the
  \code{trace.factor} are drawn in different styles so that they may
  be more easily distinguished.}
\description{}
\seealso{\code{\link{lsd.plot}}, \code{\link{margin.plot}},
  \code{\link{tukey.plot}}.}
\examples{
data(ToothGrowth)
ToothGrowth$dose = ordered(ToothGrowth$dose)
reduced <- lm(len~dose+supp, ToothGrowth)
full <- update(reduced,len~dose*supp)
margin.plot.lm(full,reduced,ToothGrowth)
}
\keyword{design}
% Converted by  version .
