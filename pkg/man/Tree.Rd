\name{Tree}
\alias{Tree}
\alias{TreeAve}
\title{Tree Weed Suppression over Time}
\usage{
data(TreeAve)
data(Tree)
}
\format{
    TreeAve data frame with 48 observations on 4 variables.
    \tabular{rlll}{
	[,1] \tab soil \tab  factor \tab   soil type \cr
	[,2] \tab reps \tab  factor \tab   replicate number \cr
	[,3] \tab score \tab numeric  \tab  germinating tomatoes of 10 \cr
	[,4] \tab height \tab numeric \tab mean plant height
    }
    Tree data frame with 96 observations on 12 variables.
    \tabular{rlll}{
	[,1] \tab soil \tab  factor \tab   soil type \cr
	[,2] \tab reps \tab  factor \tab   replicate number \cr
	[,3] \tab acc \tab factor  \tab    autoclaved or not \cr
	[,4] \tab score1 \tab numeric  \tab  germination at time 1 \cr
	[,5] \tab score2 \tab numeric  \tab  germination at time 2 \cr
	[,6] \tab score3 \tab numeric  \tab  germination at time 3 \cr
	[,4] \tab ht1 \tab numeric \tab mean plant height at time 3 \cr
	[,6] \tab score4 \tab numeric  \tab  germination at time 4 \cr
	[,4] \tab ht2 \tab numeric \tab mean plant height at time 4 \cr
	[,6] \tab score5 \tab numeric  \tab  germination at time 5 \cr
	[,4] \tab ht3 \tab numeric \tab mean plant height at time 5 \cr
	[,6] \tab score6 \tab numeric  \tab  germination at time 6
    }
}
\description{A student of plant pathology studied the `weed suppressiveness' of
tree \code{soil}s.  David Maxwell examined the sensitivity of `Rutgers'
variety tomato, known to be sensitive to biological factors, to
several different \code{soil}s.  Three of the \code{soil}s (\code{vermix,
promix} and \code{grnhous}) are in a certain sense controls, the first
two being artificial and the third coming from a controlled
environment.  A fourth \code{soil} contained only \code{grass}. The
other four are of particular interest, as the scientist wants to know
whether \code{soil} from tree \code{nm-6} is better (lower \code{score})
than those of the other three trees.  The \code{score} was the number
of germinating tomato seeds out of ten planted.  In addition, the mean
\code{height} of seedlings was measured at different times.

For each type of \code{soil}, there were six \code{soil} samples
(\code{reps(soil)}).  These \code{soil} samples were inoculated with pathogen.
Then these \code{soil} samples were split (\code{acc}), with one half
autoclaved (sterilized) and the other half not.  Thus there were
12 half-samples (\code{acc*reps(soil)}).  Each of these 12 half-samples
was placed in storage and \code{score}d on six different \code{date}s.
In addition, the mean \code{height} of seedlings was measured at three
of these times.  The primary focus is on the half-samples that were
not autoclaved.}
\source{David Maxwell, Plant Pathology Department, UW-Madison}
\references{Maxwell DL and Stanosz G (1996)
	`Assay to determine if weed suppressiveness is an inherent
	property of soil in which NM-6 is growing',
	Department of Plant Pathology, University of Wisconsin-Madison.}
\examples{
# Analysis Averaged over Time
data( TreeAve )
# soil is automatically a factor() since it has characters
TreeAve.aov <- aov( score ~ soil, TreeAve)
summary( TreeAve.aov )
#plot( TreeAve.aov )
boxplot( split( TreeAve$score, TreeAve$soil))

# Analysis over Time
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

Trees.soil <- aov(score ~ soil * date + Error(reps/soil), Trees )
# combined (sort of)
Trees.all <- aov( score ~ soil * date, Trees )

# Figure I:25.1 Tree scores over time by soil type


Tree.lsm <- lsmean( Trees.all, factors = c("soil","date") )
Tree.lsm$type <- factor( !is.na( match( Tree.lsm$soil,
   c("dtac2","nc11004","ne-332","nm-6") ) ) )

Tree.lsm <- split( Tree.lsm, Tree.lsm$type )
# following values are SE for means taken for SAS output
Tree.lsm[[1]]$se <- rep( .84086871, nrow( Tree.lsm[[1]] ) )
Tree.lsm[[2]]$se <- rep( .63433993, nrow( Tree.lsm[[2]] ) )

Tree.data <- split( Trees, Trees$type )
Tree.xlab <- c("(b) control","(a) Tree")
Tree.ylab <- c("","score means")
Tree.ylim = range( pretty( c( Tree.lsm[[1]]$pred, Tree.lsm[[2]]$pred )))

for( i in 2:1 ) { # "FALSE" comes before "TRUE"
   int.plot( Trees.all, Tree.data[[i]], c("date","soil"),
      Tree.lsm[[i]], bar.plot = "none", ylim = Tree.ylim,
      panelf = function(...) {
        se.bar( 4.5, 2, Tree.lsm[[i]]$se[1] )
        Tree.tmp <- Tree.lsm[[i]][Tree.lsm[[i]]$date==6,]
        panel.text( rep( 4.5, nrow( Tree.tmp ) ),
          0.25 + Tree.tmp$pred,
          as.character( Tree.tmp$soil ), adj = 0 )
      },
      xlab = paste( Tree.xlab[i], "soils by date" ),
      ylab = Tree.ylab[i],
      more = (i == 2), split = c(3-i,1,2,1) )
}
rm( Tree.data, Tree.xlab, Tree.ylab, Tree.lsm )

# Figure I:25.2 Tree whole plot experimental units

Tree.xlab <- c("(b) control","(a) Tree")
Tree.ylab <- c("","soil sample means")

library(lme4)
Trees.lme <- lmer(score ~ soil * date + (1|reps:soil), data = Trees )
VarCorr( Trees.lme )
anova( Trees.lme )

\dontrun{
# projections (*** CAUTION: SLOW ***)
# This works for now in Splus but not in R!
Trees.type <- aov(score ~ type/soil * date +
   Error(reps \%in\% soil \%in\% type), Trees, proj = TRUE )
Trees.proj <- proj ( Trees.type )
# whole plot
Trees.wp <- nested( Trees.proj, Trees, c("type","soil","reps"), "score",
   c(1,2), Trees.ref )
Trees.wpaov <- aov( score ~ soil, Trees.wp )
Tree.grmean <- lsmean( Trees.wpaov, fac = NULL )$pred
Tree.lsm <- lsmean( Trees.wpaov, factors = "soil" )
Tree.lsm$type <- factor( !is.na( match( Tree.lsm$soil,
   c("dtac2","nc11004","ne-332","nm-6") ) ) )
Tree.lsm <- split( Tree.lsm, Tree.lsm$type )
Tree.data <- split( Trees.wp, Trees.wp$type )
for ( i in 2:1 ) { # "FALSE" comes before "TRUE"
   Tree.lsm[[i]]$rank <- rank( -Tree.lsm[[i]]$pred )
   # have to be careful about levels of soil in type subsets
   Tree.data[[i]]$soil <- factor( Tree.data[[i]]$soil )
   attach( Tree.data[[i]] )
   names( Tree.lsm[[i]]$rank ) <- unique( soil )
   plot( jitter( Tree.lsm[[i]]$rank[soil] ), score, xaxt = "n",
      ylim = range( Trees.wp$score ),
      xlab = paste( Tree.xlab[i], "soil sample means" ),
      ylab = Tree.ylab[i] )
   detach()
   axis( 1, Tree.lsm[[i]]$rank, row.names( Tree.lsm[[i]] ) )
   abline( h = Tree.grmean, lty = 3 )
}
axis( 1, Tree.lsm[[1]]["grnhous","rank"], "grnhous" )
rm( Tree.data, Tree.xlab, Tree.ylab, Tree.grmean, Tree.lsm )

# Figure I:25.3 Tree repeated measures within soil sample

Trees.ref <- data.frame( Reference = c( Trees.proj[[1]] ))
for (i in names( Trees.proj )[-1])
   Trees.ref[[i]] <- apply( Trees.proj[[i]], 1, sum )

# split plot
Trees.sp <- nested( Trees.proj, Trees, c("type","soil","reps","date"),
   "score", c(1,3), Trees.ref )

Trees.spaov <- aov( score ~ reps*soil + date + date:soil, Trees.sp )

# split plot with time trend removed
Trees.sp$res <- resid( aov( score ~ reps*soil + date, Trees.sp ) ) +
   lsmean( Trees.spaov, fac=NULL )$pred

Trees.spres <- aov( res ~ reps*soil + date + date:soil, Trees.sp )

Tree.lsm <- lsmean( Trees.spaov, factors = c("soil","date") )
Tree.lsm$type <- factor( !is.na( match( Tree.lsm$soil,
   c("dtac2","nc11004","ne-332","nm-6") ) ) )

Tree.lsm <- split( Tree.lsm, Tree.lsm$type )
# following values are SE for contrasts taken from SAS output
Tree.lsm[[1]]$se <- rep( .42023472 / sqrt( 2 ), nrow( Tree.lsm[[1]] ) )
Tree.lsm[[2]]$se <- rep( .35818796 / sqrt( 2 ), nrow( Tree.lsm[[2]] ) )

Tree.data <- split( Trees.sp, Trees.sp$type )
Tree.xlab <- c("(b) control","(a) Tree")
Tree.ylab <- c("","pure split plot score")
for ( i in 2:1 ) { # "FALSE" comes before "TRUE"
   lsd.plot( Trees.spaov, Tree.data[[i]], c("date","soil"),
      Tree.lsm[[i]],
      xpos = 6, xlim = c(1,7.7), ylim = range( Tree.lsm$pred ),
      xaxt = "n", xlab = paste( Tree.xlab[i], "soils by date" ),
      ylab = Tree.ylab[i] )
   axis(1,1:6)
   Tree.tmp <- Tree.lsm[[i]][Tree.lsm[[i]]$date==6,]
   text( rep( 6.2, nrow( Tree.tmp ) ),
      uncollide( Tree.tmp$pred, c(.7,1.6)[i] ),
      as.character( Tree.tmp$soil ), adj = 0 )
}
rm( Tree.data, Tree.xlab, Tree.ylab, Tree.lsm, Tree.tmp )

# Figure I:25.4 Tree detrended repeated measures
Tree.lsm <- lsmean( Trees.spres, factors = c("soil","date") )
Tree.lsm$type <- factor( !is.na( match( Tree.lsm$soil,
   c("dtac2","nc11004","ne-332","nm-6") ) ) )

Tree.lsm <- split( Tree.lsm, Tree.lsm$type )
# following values are SE for contrasts take for SAS output
Tree.lsm[[1]]$se <- rep( .42023472 / sqrt( 2 ), nrow( Tree.lsm[[1]] ) )
Tree.lsm[[2]]$se <- rep( .35818796 / sqrt( 2 ), nrow( Tree.lsm[[2]] ) )

Tree.data <- split( Trees.sp, Trees.sp$type )
Tree.xlab <- c("(b) control","(a) Tree")
Tree.ylab <- c("","pure split plot score")
for ( i in 2:1 ) { # "FALSE" comes before "TRUE"
   lsd.plot( Trees.spres, Tree.data[[i]], c("date","soil"), Tree.lsm[[i]],
      xpos = 4.5, ypos = 4.25,
      xlim = c(1,7.7), ylim = range( Tree.lsm$pred ),
      xaxt = "n",
      xlab = paste( Tree.xlab[i], "soils by date" ),
      ylab = Tree.ylab[i] )
   axis(1,1:6)
   Tree.tmp <- Tree.lsm[[i]][ Tree.lsm[[i]]$date == 6, ]
   text( rep( 6.2, nrow( Tree.tmp ) ),
      uncollide( Tree.tmp$pred, c(.3,.6)[i] ),
      as.character( Tree.tmp$soil ), adj = 0 )
}
rm( Tree.data, Tree.xlab, Tree.ylab, Tree.lsm, Tree.tmp )

# Figure I:26.1 Tree partial correlations over time
Tree.tmpc <- length( unique( Trees.sp$date ))
Tree.tmpr <- nrow( Trees.sp ) / Tree.tmpc
Tree.tmp <- matrix( resid( Trees.spaov ), Tree.tmpr, Tree.tmpc )
Tree.tmpid <- c("C","T")[ codes( Trees.sp$type[1:Tree.tmpr] ) ]

splom( Tree.tmp, varnames = as.character( unique( Trees.sp$date )),
   groups = Tree.tmpid )
rm( Tree.tmp, Tree.tmpr, Tree.tmpc )
}
} 
% Converted by  version .
\keyword{datasets}
