\name{pdalme}
\alias{pdalme}
\title{PDA LME Code Checking}
\usage{
}
\description{
Trying to figure out how best to use library lme4
}
\examples{
library( lme4 )
########################################################3
data( Berry )
Berry.lme <- lmer( conc ~ tissue * date + (1|sample),
   data = Berry, method = "REML" )
print( Berry.lme )
anova( Berry.lme )
########################################################3
data( Brassica )
# drop Crop 2
Brassica <- data.frame(Brassica[ Brassica$crop != 2, ])
Brassica$crop <- factor( as.character( Brassica$crop ))
Brassica$yrblk = interaction(Brassica[,c("year","block")])
Brassica$cropyrblk = interaction(Brassica[,c("crop","year","block")])

Brassica1 <- data.frame( Brassica[ Brassica$year == 1, ] )
Brassica.lme <- lmer(yldkga ~ crop * pyrrat + (1|block:crop),
   data = Brassica1 )
anova( Brassica.lme )

Brassica.fulle = lmer(yldkga ~ crop * pyrrat * year +
  (1|yrblk:cropyrblk), data = Brassica )
anova(Brassica.fulle)
Brassica.rede = lmer(yldkga ~ crop * year + pyrrat + crop:pyrrat +
  (1|yrblk:cropyrblk), data = Brassica )
anova(Brassica.rede)
anova(Brassica.rede,Brassica.fulle)
## plot does not work with lmer yet
lsd.plot(Brassica.rede,factor=c("pyrrat","crop"),
   xlab = "pyr rate by crop",
   ylab = "mean crop yield" )
########################################################
data( Drink )
Drink2 <- Drink[Drink$period < 3,]
Drink4 <- Drink[!is.na(match(Drink$cow,Drink$cow[Drink$period == 4])),]

Drink$cow <- factor( Drink$cow )
Drink$month <- ordered( Drink$month )
Drink$period <- ordered( Drink$period )

Drink2$cow <- factor( Drink2$cow )
Drink2$month <- ordered( Drink2$month )
Drink2$period <- ordered( Drink2$period )

Drink4$cow <- factor( Drink4$cow )
Drink4$month <- ordered( Drink4$month )
Drink4$period <- ordered( Drink4$period )

Drink2.lme <- lmer(milk ~ period * water + (1|cow),
  data = Drink2 )
summary(Drink2.lme)
anova(Drink2.lme)
VarCorr( Drink2.lme)

Drink4.lme <- lmer(milk ~ period * water + (1|cow),
  data = Drink4 )
Drink4.lme
anova(Drink4.lme)
VarCorr( Drink4.lme)

int.plot( Drink2.lme, factors = c("period","water"),
   xpos = 1.25, xlab = "(a) cows in first two periods",
   ylab = "7-day milk yield (lb)", bar = "ellipse" )

int.plot( Drink4.lme, factors = c("period","water"),
   xpos = 1.75, xlab = "(b) cows in all four periods",
   ylab = "", bar = "ellipse" )
int.plot( Drink4.lme, factors = c("water","period"),
   xpos = 1.75, xlab = "(c) cold vs hot",
   ylab = "", bar = "ellipse" )
########################################################
data( Random )
Random.frame <- data.frame( response = c( as.matrix( Random )))
Random.frame$class <- factor( rep( 1:10, rep( 10,10 )))
Random.lme <- lmer(response ~ 1 + (1|class), Random.frame )
summary( Random.lme )
VarCorr( Random.lme )
########################################################
data( Rantwo )
Rantwo$seed = factor( Rantwo$seed )
Rantwo$farm = factor( Rantwo$farm )
Rantwo$seedfarm = with( Rantwo, seed:farm)
Rantwo.lme <- lmer(yield ~ 1 + (1|seed) + (1|farm) + (1|seedfarm),
   data = Rantwo )
summary( Rantwo.lme )
Rantwo.lmep <- lmer(yield ~ 1 + (1|seed) + (1|farm),
   data = Rantwo )
anova(Rantwo.lme,Rantwo.lmep)
########################################################
data( Tree )
Tree <- Tree[ Tree$acc == "notauto", ]
Tree$acc <- NULL
Tree$type <- factor( !is.na( match( Tree$soil,
   c("dtac2","nc11004","ne-332","nm-6") )))
Tree$reps <- factor( Tree$reps )

Trees <- Tree[ rep( seq( nrow( Tree ) ), 6 ), c("soil","reps","type") ]
Trees$date <- ordered( rep( 1:6, rep( nrow( Tree ), 6 ) ) )
Trees$score <- unlist( Tree[ , paste( "score", 1:6, sep="" ) ] )
row.names(Trees) <- interaction( Trees[, c("soil","reps","date") ] )
Trees.lme <- lmer(score ~ soil * date + (1 | reps:soil ), data = Trees )
anova(Trees.lme)
VarCorr(Trees.lme)
########################################################
data( WheatPDA )
class( WheatPDA$ploidy ) = factor( WheatPDA$ploidy, ordered = FALSE )
WheatPDA$sp = interaction( WheatPDA$ploidy, WheatPDA$species, drop = TRUE)
WheatPDA$cr = interaction( WheatPDA$sp, WheatPDA$cross, drop = TRUE)
WheatPDA$ac = interaction( WheatPDA$cr, WheatPDA$accession, drop = TRUE)
\dontrun{
## following lead to singular X'X matrices
Wheat.lme = lmer( tr1 ~ ploidy + sp + (1|cr) + (1|ac),
  data = WheatPDA, na.action = na.omit )
Wheat.lme = lmer( tr1 ~ ploidy + sp + (1|cross),
  data = WheatPDA, na.action = na.omit )
}
}
% Converted by  version .
\keyword{datasets}
