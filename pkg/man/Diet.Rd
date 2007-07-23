\name{Diet}
\alias{Diet}
\title{Diet Effect on Dry Matter Intake}
\usage{
data(Diet)
}
\format{
    Diet data frame with 12 observations on 7 variables.
    \tabular{rlll}{
	[,1] \tab trt \tab factor  \tab treatment identifier \cr
	[,2] \tab bwt \tab numeric \tab body weight \cr
	[,3] \tab bwg \tab numeric \tab body weight gain \cr
	[,4] \tab fi  \tab numeric \tab diet intake \cr
	[,5] \tab fe  \tab numeric \tab diet(?) \cr
	[,6] \tab cla \tab numeric \tab CLA amount \cr
 	[,6] \tab la  \tab numeric \tab LA amount
    }
}
\description{This experiment (Wattiaux, Combs and Shaver 1994) involved 60 cows
initially, although some were lost 
during the study.  The cows were randomly assigned to one of 6 diets
and followed for a number of weeks.  Diets were begun after the third
week, allowing the animals some initial time to adjust to their new
environment.  Interest focuses on the effect of diet on the average
dry matter intake (\code{dmi}), the amount of food eaten by each cow.
The data you have are a baseline value (\code{covar} = dmi
for week 3), average \code{dmi} during subsequent weeks, the
number of subsequent \code{weeks}.  In addition there are cow and
block identifiers.

The experiment was a randomized block design.  That is, cows were \code{blocked} by time, the first 6 cows were randomly assigned among the 6 diets,
and so on.  ``Proper'' analysis should take account of this, along with
the initial capacity of each animal (its \code{covariate} of dry matter
intake at 3 weeks) and possibly weighting by the number of \code{weeks}
each cow was on trial.}
\source{Michel Wattiaux,
	International Agriculture Program, CALS, UW-Madison}
\references{Wattiaux MA, Combs DK and Shaver RD (1994)
	``Lactational responses to ruminally undegradable protein by
	dairy cows fed diets based on alfalfa silage'',
	\emph{J Dairy Science 77}, 1604-1617.}
\examples{
data( Diet )

Diet$trt <- factor(Diet$trt)
Diet$block <- factor(Diet$block)
#
# ordinary anova
Diet.fit <- aov(dmi ~ trt, Diet)

# ancova with covar = first 3 weeks dmi
Diet.cov <- lm(dmi ~ covar + trt, Diet, na.action=na.omit)

# ancova with blocking and weighting by number of weeks on trial
Diet.blk <- lm(dmi ~ covar + block + trt, Diet,
   weights=Diet$weeks, na.action=na.omit)

# Figure F17.2: Diet Adjusted for Initial Capacity
## regression lines in following plot are not parallel
Diet = Diet[!is.na(Diet$covar), ]
Diet$pred = predict(Diet.cov)
trellis.par.set(theme=col.whitebg())
print(xyplot(dmi ~ covar, Diet, groups = trt, type = "p",
  pch = levels(Diet$trt), cex = 1.5,
  col = 1 + as.numeric(levels(Diet$trt)),
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    pred = predict(Diet.cov)
    for(i in levels(Diet$trt)) {
      ii = Diet$trt == i
      i = as.numeric(i)
      panel.lines(Diet$covar[ii], Diet$pred[ii],
        lty = i, col = i + 1)
    }
    ## SD bar
    lsd.bar(Diet.cov, xpos = 24, ypos=20, cap = "SD")
  },
  xlab = "covariate = first 3 weeks",
  ylab = "dry matter intake (dmi)",
  main = "Figure F:17.2 Diet Adjusted for Initial Capacity"))
data.frame(se.bar(26, 20, std.dev(Diet.cov), cap = "SD"))
} 
% Converted by  version .
\keyword{datasets}
