\name{trans.plot}
\alias{trans.plot}
\title{Diagnostic plot to find power transform to remove
    interaction. See Emerson chapter in Hoaglin, Mosteller and Tukey
    (1991, ch. 13E) Fundamentals of Exploratory Analysis of Variance
    (Wiley);or see Mosteller and Tukey (1977, ch. 9E,F) Exploratory Data
    Analysis (Addison-Wesley).}
  \usage{
trans.plot(fit, predictors=all.predictors(fit)[1:2],
  terms=predict(fit, type = "terms"), xlab, ylab, \dots)
}
\description{}
\arguments{
  \item{fit}{linear model (lm) fit}
  \item{predictors}{predictors to compare.  They must be in the model, as
    well as theirinteraction in the specified order.}
  \item{terms}{estimates of model terms}
  \item{xlab,ylab}{x and y axis labels}
  \item{\dots}{Plot parameters}
}
\value{Power = 1 - slope estimate.  SE comes from simple linear
    regression and is at best advisory.  Consider tranformation $y^p$ =
    exp(p*log(y)), or transform log(y) if p=0.  In practice, examine
    nearest half integer.}
\seealso{\code{\link{effect.plot}}.}
\examples{
\dontrun{
# Bekk data from Hoaglin, Mosteller and Tukey (1991)
bekk.fit <- lm(smooth~lab*mat, bekk)
effect.plot(bekk.fit)trans.plot(bekk.fit)
bekk.lfit <- lm(log(smooth)~lab*mat, bekk)
effect.plot(bekk.lfit)
trans.plot(bekk.lfit)
}
}
\keyword{hplot}
% Converted by  version .
