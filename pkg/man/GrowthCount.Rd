\name{GrowthCount}
\alias{GrowthCount}
\title{Counting Data for Growth Medium}
\usage{
data(Growth)
}
\description{See Details for \code{\link{GrowthCount}}.}
\source{Sue Stieve \& Dennis Stimart}
\references{Venables and Ripley (1999) \emph{Modern Applied Statistics
	with S-PLUS}, 3rd ed., ch. 7.}
\seealso{\code{\link{Growth}}.}
\examples{
data( Growth )

Growth <- Growth[,c("trt","advplt")]
# drop the added BHTA control ( trt = 20 )
Growth <- Growth[ Growth$trt != 20, ]
 
Growth$BA <- c(rep(0,4),rep(0.44,4),rep(4.4,3),rep(NA,9),0)[1+Growth$trt]
Growth$TDZ <- c(rep(c(0,0.2,2,20),3)[1:11],rep(NA,9),0)[1+Growth$trt]
Growth$code <- c(0:9,"a",rep(NA,9),"C")[1+Growth$trt]

Growth$trt <- factor(Growth$trt)
Growth$BA <- factor(Growth$BA)
Growth$TDZ <- factor(Growth$TDZ)
Growth <- Growth[!is.na(Growth$advplt),]

sample.size( Growth$BA, Growth$TDZ )
hist( Growth$advplt, nclass = 25 )
table( Growth$advplt )
Growth$nonzero <- as.numeric( Growth$advplt > 0 )

Growth.bin <- glm( nonzero ~ BA * TDZ, data = Growth,
   family = binomial )
anova( Growth.bin, test = "Chisq" )
Growth.bin <- glm( nonzero ~ TDZ *  BA, data = Growth,
   family = binomial )
anova( Growth.bin, test = "Chisq" )
drop1( Growth.bin, formula( Growth.bin), test = "Chisq" )

Growth.des <- data.frame(
   BA = factor( rep( levels( Growth$BA ), 4 )),
   TDZ = factor( rep( levels( Growth$TDZ ), rep(3,4) )))
Growth.pred <- 1 / ( 1 + exp( - predict(
   Growth.bin, Growth.des )))
plot( c(1,21),range( Growth.pred ), type = "n", log = "x",
   xlab = "TDZ", ylab = "fraction nonzero" )
for( i in levels( Growth.des$BA ))
{
   Growth.BA <- i == Growth.des$BA
   text( 1 + unfactor( Growth.des$TDZ[Growth.BA] ),
      Growth.pred[Growth.BA], i )
}
rm( Growth.BA )

Growth.pois <- glm( advplt ~ TDZ * BA, data = Growth,
   family = poisson )
anova( Growth.pois )
Growth.apois <- glm( advplt ~ TDZ + BA, data = Growth,
   family = poisson )
anova( Growth.apois )
drop1( Growth.apois, formula( Growth.apois ), test = "F" )
} 
% Converted by  version .
\keyword{datasets}
