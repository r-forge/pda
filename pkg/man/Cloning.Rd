\name{Cloning}
\alias{Cloning}
\alias{ClonCrit}
\title{Cloning virus in plant sap}
\usage{
data(Cloning)
data(ClonCrit)}
\format{
    Cloning data frame with 64 observations on 5 variables.
    \tabular{rlll}{
	[,1] \tab clone \tab   factor \tab  clone identifier\cr
	[,2] \tab reps \tab      factor \tab  replication number \cr
	[,3] \tab titer \tab  numeric \tab concentration of virus\cr
	[,4] \tab code \tab  factor \tab plot code\cr
	[,5] \tab type \tab   factor \tab cell type
    }
    ClonCrit data frame with 20 observations on 6 variables.
    \tabular{rlll}{
	[,1] \tab set \tab   factor \tab  all data or balanced \cr
	[,2] \tab means \tab   numeric \tab number of means \cr
	[,3] \tab SNK \tab      numeric \tab  Student-Neuman-Kuels \cr
	[,4] \tab DUNCAN \tab  numeric \tab Duncan's Multiple Range\cr
	[,5] \tab REGWQ \tab  numeric \tab REGWQ\cr
	[,6] \tab REGWF \tab   numeric \tab REGWFtype
    }
}
\description{A plant scientist measured the concentration of a particular virus in
plant sap using ELISA (enzyme-linked immunosorbent assay)
(Novy 1992).  The study included \code{13} potato \code{clones}
(plants reproduced asexually), with \code{2} commercial
cultivars, \code{5} somatic hybrids, \code{5} progeny of the somatic
hybrids, and one clone of \emph{Solanum etuberosum} (a species related
to potato). Of the \code{5} progeny of the somatic hybrids, two were
classified as susceptible and three as resistant to the virus.   The
scientist wants to understand the resistance to the virus among
these \code{13 clones}.  Plant sap was taken from \code{5} inoculated
plants of each \code{clone},  for a total of \code{65} measurements of
titer.  Unfortunately, one measurement was lost during processing of
the samples.  Figure <A HREF="clondat.s">clondat.s</A>
shows the raw data by \code{clone}.
The (a) figure, ordered by \code{clone} identifier, indicates considerable
range in center and spread but is otherwise not very informative.  The
(b) figure orders by \code{clone} mean, jittered to
mitigate overlap, using a separate plot symbol for each
\code{clone}. The dashed identity line helps track the spread around the
mean.  This latter plot suggests there is less spread in titer for
\code{clones} at the extremes than for those in the middle of the range.
}
\source{Professor Richard G Novy (\url{mailto:novy@prairie.nodak.edu}),
	Plant Sciences Department, ND St U}
\references{Novy RG (1992) ``Characterization of somatic hybrids between
\emph{Solanum etuberosum} and diploid, tuber-bearing Solanums.
PhD Dissertation, Department of Plant Pathology, UW-Madison. 123 pp.}
\examples{
data( Cloning )

Cloning <- Cloning[ !is.na(Cloning$titer), ]
Cloning$clone <- factor(Cloning$code)

# ANOVA table
Cloning.fit <- aov( titer ~ clone, Cloning )
summary( Cloning.fit )

# estimate means
attach( Cloning )
Cloning.means <- data.frame( clone=levels(code),
   pred=c( tapply( titer, code, mean )))
Cloning.means$sd <- c( tapply( titer, code,
   function( x ){ sqrt( var( x ) ) } ))
Cloning.means$n <- c( tapply( titer, code, length ))
Cloning.means$se <- Cloning.means$sd /
   sqrt( Cloning.means$n )
Cloning.means$type <- type[ !duplicated( code ) ]
detach( )

# or try the PDA routine:
Cloning.lsm <- lsmean( Cloning.fit )

# estimate common SD
std.dev( Cloning.fit )

# preselected contrasts
aov(titer ~ clone, Cloning, 
   subset = is.na( match( Cloning$type, c("cult","etb","odd") )))

# critical values for multiple comparisons (from SAS output)
Cloning.crit <- data.frame( all=c(306,555,737,530,278),
   bal=c(254,435,527,414,232),
   row.names=c("LSD","BON","SCHEFFE","HSD","WALLER") )
for (i in names(Cloning.crit))
   names(Cloning.crit[[i]]) <- row.names(Cloning.crit)

# critical values for multi-stage testing (MST) multiple comparisons
data( ClonCrit )
Cloning.tmp <- ClonCrit$set
ClonCrit$set <- NULL
ClonCrit <- split(ClonCrit,Cloning.tmp)
rm(Cloning.tmp)

# Figure B:4.3 Cloning Data by Group
print(stripplot(titer ~ clone, Cloning, jitter = TRUE,
  panel = function(x,y,...) {
    panel.abline(v=unclass(x),lty=3,col="grey")
    panel.stripplot(x,y,...)
  },
  xlab="(a) clone identifier",ylab="titer",
  main = "Figure B:4.3" ),
  more = TRUE, split = c(1,1,2,1))

Cloning$pred <- jitter( predict( Cloning.fit))
print( xyplot( titer ~ pred, Cloning, groups = code,
  type = "p", cex = 1, pch = levels( Cloning$code ),
  panel = function(x,y,...) {
    panel.superpose(x,y,...)
    panel.abline( 0, 1, lty = 2 )
  },
  xlab="(b) clone mean",ylab="titer",
  main = "Cloning Data by Group" ),
  split = c(2,1,2,1))

# Figure B:4.4 Plot of Clone Mean vs SD
print( xyplot( sd ~ pred, Cloning.means, groups = clone,
  type = "p", pch = levels( Cloning.means$clone ), cex = 2,
  xlab="(a) mean by clone identifier", ylab="clone SD",
  main = "Figure B:4.4" ),
  more = TRUE, split = c(1,1,2,1))
print( xyplot( sd ~ pred, Cloning.means, groups = type,
  type = "p", pch = levels( Cloning.means$type ), cex = 2,
  xlab="(b) mean by type", ylab="clone SD",
  main = "Clone Mean vs SD" ),
  split = c(2,1,2,1))

# Figure B:4.7 95\% Critical Values by Groups
ci.plot( Cloning.fit, lsm = Cloning.lsm,
  rdf=df.resid( Cloning.fit), sort=TRUE, ylim = c(-150,2000),
  xlab = "", ylab = "(a) with common S",
  main = "Figure B:4.7",
  panelf = function(...) {
    panel.abline( h = 0, lty = 2 )
  },
  more = TRUE, split = c(1,1,2,1) )

ci.plot( Cloning.fit, lsm = Cloning.means,
  rdf = Cloning.means$n-1, sort=TRUE, ylim = c(-150,2000),
  xlab = "", ylab = "(b) with group SDs",
  main = "Cloning 95 \% CIs",
  panelf = function(...) {
    panel.abline( h = 0, lty = 2 )
  },
  split = c(2,1,2,1) )

# Figure B:6.2 Cloning Unbalanced Multiple Comparisons
# critical values for Cloning data with all 13 groups
Cloning.x <- c( 2, max( ClonCrit$all$means ))
Cloning.y <- range( unlist( c( Cloning.crit$all,
   ClonCrit$all[,-1] )))
Cloning.tmp <- round( mean( Cloning.x ) )
Cloning.zaxis <- seq( 1.5, 5, .5 )
Cloning.z <- Cloning.zaxis * sqrt(2/5)*239

plot( Cloning.x, Cloning.y, type="n",
   xlab="(a) based on F test", ylab="" )
title( "Figure B:6.2(a) Unbalanced Comparisons" )
axis( 4, Cloning.z,
   labels = as.character( Cloning.zaxis ))
for ( i in c("LSD","BON","SCHEFFE","WALLER") ) {
	abline( Cloning.crit[i,"all"], 0, lty=2 )
	text( Cloning.tmp, Cloning.crit[i,"all"], i )
}
for ( i in c("REGWF") ) {
	lines( ClonCrit$all$means, ClonCrit$all[[i]] )
	text( Cloning.tmp, ClonCrit$all[[i]][Cloning.tmp], i )
}

plot( Cloning.x, Cloning.y, type="n",
   xlab = "(b) based on range", ylab = "critical value" )
title( "Figure B:6.2(b) Unbalanced Comparisons" )
axis( 4, Cloning.z, labels = as.character( Cloning.zaxis ))

Cloning.p <- NULL
for ( i in c("LSD","HSD") )
	Cloning.p[i] <- Cloning.crit[i,"all"]
for ( i in c("REGWQ","SNK","DUNCAN") )
	Cloning.p[i] <- ClonCrit$all[[i]][Cloning.tmp]
Cloning.p["SNK"] <- Cloning.p["SNK"] - 25
Cloning.p["HSD"] <- Cloning.p["HSD"] + 25
for ( i in c("LSD","HSD")) {
	abline( Cloning.crit[i,"all"], 0, lty=2 )
	text( Cloning.tmp, Cloning.p[i], i )
}
for ( i in c("REGWQ","SNK","DUNCAN") ) {
	lines( ClonCrit$all$means, ClonCrit$all[[i]] )
	text( Cloning.tmp, Cloning.p[i], i )
}
rm( Cloning.tmp, Cloning.x, Cloning.y, Cloning.p, Cloning.z )

# Figure B: Critical Values for Comparison of 9 Clones
# critical values for Cloning data with 9 groups
Cloning.x <- c( 2, max( ClonCrit$bal$means) )
Cloning.y <- range( unlist( c( Cloning.crit$bal,
   ClonCrit$bal[,-1] ) ) )
Cloning.tmp <- round( mean( Cloning.x ) )
Cloning.z <- Cloning.zaxis * sqrt( 2 / 4.91) * 198

plot( Cloning.x, Cloning.y, type="n",
   xlab="(a) based on F test", ylab="" )
title( "Figure B:6.3(a) Balanced Comparisons" )
axis(4, Cloning.z, lab = as.character( Cloning.zaxis ) )
for ( i in c("LSD","BON","SCHEFFE","WALLER") ) {
	abline( Cloning.crit[i,"bal"], 0, lty=2 )
	text( Cloning.tmp, Cloning.crit[i,"bal"], i )
}
for ( i in c("REGWF") ) {
	lines( ClonCrit$bal$means, ClonCrit$bal[[i]] )
	text( Cloning.tmp, ClonCrit$bal[[i]][Cloning.tmp], i )
}

plot( Cloning.x, Cloning.y, type="n",
   xlab="(b) based on range", ylab="critical value" )
title( "Figure B:6.3(b) Balanced Comparisons" )
axis(4, Cloning.z, labels = as.character( Cloning.zaxis ))
Cloning.p <- NULL
for ( i in c("LSD","HSD") )
	Cloning.p[i] <- Cloning.crit$bal[i]
for ( i in c("REGWQ","SNK","DUNCAN") )
	Cloning.p[i] <- ClonCrit$bal[[i]][Cloning.tmp]
Cloning.p["SNK"] <- Cloning.p["SNK"] - 25
Cloning.p["HSD"] <- Cloning.p["HSD"] + 25
for ( i in c("LSD","HSD") ) {
	abline( Cloning.crit[i,"bal"], 0, lty=2 )
	text( Cloning.tmp, Cloning.p[i], i )
}
for ( i in c("REGWQ","SNK","DUNCAN") ) {
	lines( ClonCrit$bal$means, ClonCrit$bal[[i]] )
	text( Cloning.tmp, Cloning.p[i], i )
}
rm( Cloning.tmp, Cloning.x, Cloning.y, Cloning.p, Cloning.z,
   Cloning.zaxis )

# Figure B:6.4 Critical Value by Number of Groups

Cloning.groups <- 2:20
Cloning.alpha <- .05

plot( range(Cloning.groups), c(1,6), type="n",
	xlab="number of groups", ylab="critical value" )
title( "Figure B:6.4 Critical Value by Number of Groups" )
axis(1,c(2,13))
# SCHEFFE
lines( Cloning.groups, sqrt( qchisq( 1 - Cloning.alpha,
   Cloning.groups - 1 ) ), lty=2 )
text( 13, 5, "SCHEFFE" )
# TUKEY's HSD
lines( Cloning.groups,
   c(2.77,3.32,3.63,3.86,4.03,4.17,4.29,4.39,4.47,4.55,4.62,
   4.68,4.74,4.8,4.84,4.89,4.93,4.97,5.01) / sqrt(2), lty=3 )
text( 13, 3, "HSD" )
# BON
Cloning.comp <- Cloning.groups * ( Cloning.groups - 1 )
lines( Cloning.groups,
   qnorm( 1 - Cloning.alpha / Cloning.comp ) )
text( 13, 3.8, "BONFERRONI" )
# LSD
lines( c(2,20), rep( qnorm( 1 - Cloning.alpha / 2), 2 ) )
text( 13, 1.6, "LSD" )

rm( Cloning.alpha, Cloning.groups, Cloning.comp )

# Figure B:6.5 Bonferroni Comparisons to Match Scheffe' Method

Cloning.groups <- 2:20
Cloning.alpha <- .05

# equivalent comparisons between B and S
Cloning.a <- log( Cloning.alpha / ( 2 * ( 1 - pnorm( sqrt(
   qchisq( .95, Cloning.groups - 1 ))))))
plot(Cloning.groups, Cloning.a, type="l", yaxt="n",
	xlab="number of groups", ylab="comparisons" )
title( "Figure B:6.5 Match of Bonferroni to Scheffe" )
axis( 1, c(2,13) )
axis( 2, log( c(1,10,78,1000,10000,100000,1000000) ),
	c("1","10","78","1000","10^4","10^5","10^6"), adj=1 )
axis( 2, log(100), FALSE )
lines( c(13,13), c(0,Cloning.a[12]), lty=2 )
lines( c(2,13), rep(log(78),2), lty=2 )

rm( Cloning.groups, Cloning.alpha, Cloning.a )
}
% Converted by  version .
\keyword{datasets}
