\name{lsmean}
\alias{lsmean}
\alias{lsmean.default}
\alias{lsmean.lm}
\alias{lsmean.lme}
\alias{lsmean.lmer}
\alias{lsmean.listof}
\title{Least Squares Means}
\usage{
lsmean(object, \dots)
\method{lsmean}{default}(object, \dots, factors, effects = FALSE, se.fit = TRUE,
   adjust.covar = TRUE)
\method{lsmean}{lm}(object, data, factors, expr, contrast, effects = FALSE,
   se.fit = TRUE, adjust.covar = TRUE, pdiff = FALSE, reorder = FALSE,
   lsd, level = .05, rdf, coef, cov, \dots)
\method{lsmean}{lme}(object, data, factors, \dots, rdf, coef, cov)
\method{lsmean}{lmer}(object, data, factors, expr, \dots, rdf, coef, cov)
\method{lsmean}{listof}(object, data, factors, stratum, expr, contrast, \dots)
}
\arguments{
  \item{object}{response vector (default) or model object (lm).}
  \item{\dots}{factors and covariates (must be same length as y).}
  \item{data}{data frame in which to interpret variables(found from
    object if missing).}
  \item{factors}{character vector containing names of x.factor and
    trace.factoras first two entries.  Must be in \code{names(data)} and
    \code{labels(object)}.Default is all factor names.}
  \item{effects}{drop intercept if \code{TRUE} (only works properly with
    sum-to-zero contrasts).}
  \item{se.fit}{compute pointwise standard errors if T.}
  \item{adjust.covar}{adjust means to average covariate values if
    T; otherwise use covariate mean for each combination of
    factors.}
  \item{pdiff}{Include letters to signify significant differences.}
  \item{reorder}{Reorder means from largest to smallest.}
  \item{lsd}{Include average LSD if \code{TRUE} (also need \code{pdiff=TRUE}).}
  \item{level}{Significance level for \code{pdiff} calculations.}
  \item{rdf}{Residual degrees of freedom.}
  \item{coef}{Coefficients for fixed effects in object.}
  \item{cov}{Covariance matrix for fixed effects.}
  \item{expr}{Call expression (formula)}
  \item{contrast}{Type of contrasts (default is attribute
    \code{contrasts} of \code{object})}
  \item{stratum}{Name of stratum for lsmean calculation as character string.}
}
\description{Caution: This routine is not fully tested for models with nested
  factors or mixed models. Please check results against another
  package (e.g. SAS proc mixed). It appears to correctly handle \code{lme}
  objects, but does
  not work well for \code{aov} objects that include \code{Error()}
  type nesting in the formula. Further, it does not properly handle
  polynomial terms--only the linear term is included. For now, create
  dummies like x2 = x*x manually and include x2 in your model. 
}
\value{Data frame containing unique factor levels of factors, predicted
  response (pred) and standard errors (se). WARNING: lsmean may not
  function properly if there are empty cells. Standard errors for mixed
  models using methods \code{lmer} and \code{listof} are not fully
  debugged.
}
\seealso{\code{\link[base]{predict}}.}
\examples{
\dontrun{
lsmean(y,x1,x2)
# the following does the same thing
fit <- lm(y~x1+x2)
data <- data.frame(y,x1,x2)
lsmean(fit,data,factors=c("x1","x2")
}
}
\keyword{design}
% Converted by  version .
