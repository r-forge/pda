\name{box.fig}
\alias{box.fig}
\title{draw boxplot figure at specified position}
\usage{box.fig(fit, xpos=1, width=0.2, jitter=width)}
\arguments{
    \item{fit}{linear model (lm) fit}
    \item{xpos}{horizontal (x) position of boxplot}
    \item{width}{width of boxplot}
    \item{jitter}{uniform jitter of outliers, if any}}
\value{list as returned from boxplot(fit,plot=F)}
\description{}
\seealso{\code{\link[base]{boxplot}},
    \code{\link{effect.plot}}.}
\examples{
\dontrun{
pred <- fitted(fit)
plot(pred,resid(fit),xlim=range(pred)*c(1,1.1))
box.fig(fit,max(pred)*1.1)
}
}
\keyword{hplot}
% Converted by  version .
