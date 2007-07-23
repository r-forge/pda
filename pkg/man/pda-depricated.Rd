\name{pda-depricated}
\alias{pda-depricated}
\alias{char}
\alias{choose}
\alias{cor.coef}
\alias{cov.coef}
\alias{dotna}
\alias{fact}
\alias{get.list}
\alias{ghostview}
\alias{harmonic.mean}
\alias{mlines}
\alias{mplot}
\alias{mpoints}
\alias{mqqnorm}
\alias{perms}
\alias{plotter}
\alias{se.coef}
\alias{split.plot}
\title{Depricated pda functions.}
\description{
  These were old and in the way.
}
\usage{
char()
choose()
cor.coef()
cov.coef()
dotna()
fact()
get.list()
ghostview()
harmonic.mean()
mlines()
mplot()
mpoints()
mqqnorm()
perms()
plotter()
se.coef()
split.plot()
}
\details{
  \item{char}{Use \code{pch} in plot routines.}
  \item{cor.coef,cov.coef,se.coef}{Use \code{\link[base]{vcov}} or \code{\link[base]{VarCorr}}.}
  \item{dotna}{Use \code{na.strings} in \code{read.table}.}
  \item{choose,fact,perms}{See \code{\link[base]{Special}} Functions of Mathematics.}
  \item{get.list}{Use \code{data[factors]}.}
  \item{ghostview}{Use system call if needed.}
  \item{harmonic.mean}{Use \code{1/mean(1/data)}.}
  \item{mlines,mplot,mpoints}{Use \code{\link[lattice]{xyplot}}.}
  \item{mqqnorm}{Use \code{\link[lattice]{qqmath}}.}
  \item{plotter}{Use \code{help(\link[grdevices]{Devices})} to find available devices.}
  \item{split.plot}{Use \code{\link{splitplot}}.}
}
\author{ Brian S. Yandell, yandell@stat.wisc.edu }
\keyword{internal}
