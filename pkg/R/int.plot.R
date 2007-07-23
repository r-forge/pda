#####################################################################
##
## $Id: int.plot.R,v 1.1 2005/03/24 yandell@stat.wisc.edu Exp $
##
##     Copyright (C) 2005 Brian S. Yandell
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the
## Free Software Foundation; either version 2, or (at your option) any
## later version.
##
## These functions are distributed in the hope that they will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## The text of the GNU General Public License, version 2, is available
## as http://www.gnu.org/copyleft or by writing to the Free Software
## Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
##
##################################################################
### to do: 24 mar 2005
### remove mlines and mplot everywhere (see bwplot in lattice)
###   check man/*.Rd files
###   troubleshoot int.plot.lm
### fix kludges in lsmean
##################################################################
lsd.plot <- function(object, ...) UseMethod("int.plot")
int.plot <- function(object, ...) UseMethod("int.plot")

int.plot.default <- function(object, response, group, 
                             xlab=xlabs, ylab=ylabs, ...)
{
  xlabs <- deparse(substitute(object))
  ylabs <- deparse(substitute(response))

  if (!is.factor(object))
    object <- factor(object)
  if( missing(group) ) {
    data <- data.frame(response, object)
    names(data) <- c(ylabs, xlabs)
    formula.y <- formula( paste( ylabs, xlabs, sep="~" ))
    fit <- lm(formula.y, data)
  }
  else {
    glabs <- deparse(substitute(group))
    if (!is.factor(group))
      group <- factor(group)
    data <- data.frame(response, object, group)
    names(data) <- c(ylabs, xlabs, glabs)
    formula.y <- formula( paste( ylabs, ".^2", sep="~" ))
    fit <- lm(formula.y, data)
  }
  invisible(int.plot.lm(fit, data, xlab=xlab, ylab=ylab, ...))
}
int.plot.lmer <- int.plot.lme <-
  function( object, data = eval( object@call$data ),
           ... ) int.plot.lm( object, data, ... )

## NB: lsd and pdiff need only be TRUE for some settings in lsmean call.
int.plot.lm <- function(object, data=eval(object$call$data), 
                        factors=all.factors(object, data), 
                        lsm=lsmean(object, data, factors, pdiff=TRUE, lsd=TRUE), 
                        ci=ci.width(object, data, factors, lsm, rdf, level), 
                        offset=0.5, bar.plot="lsd", 
                        ylim=ylims, type="b", 
                        width=.01, sort.mean=FALSE, edge=wi, 
                        xlab=xlabs, ylab=as.character(formula(object)[[2]]), 
                        level=.05, rdf = df.resid(object), 
                        fine=.01, xpos=xposn, ypos=yposn, cex=2, white=TRUE,
                        panelf = NULL, ...)
{
  myaxis = list(...)$xyplot
  is.xyplot = !is.null( myaxis )
  if( !is.xyplot )
    myaxis = function(...){}
  if( is.null( panelf ))
    panelf = function(...){}
  
  bar.plot <- c("lsd","ci","test","ellipse","none")[
                pmatch( bar.plot, c("lsds","cis","tests","ellipses"),
                nomatch = 5 ) ]

  ## make room for confidence intervals or test intervals
  if( pmatch( bar.plot, c("ci", "test"), nomatch = 0 )) {
    wi <- .02
    if ( bar.plot == "test" )
      ci <- ci / sqrt(2)
    ylims <- c(min(lsm$pred-ci), max(lsm$pred+ci))
  }
  else {
    wi <- 0
    ylims <- range(lsm$pred)
  }
  
  ## set up plotting groups if any
  xlabs <- factors[1]
  lf <- length(factors)
  if (lf == 1) { # only 1 factor
    group <- NULL
    x.factor <- lsm[[factors]]
  }
  else {
    if (lf > 2) { # x is 1st factor, all others have traces
      group <- as.character(interaction( lsm[, factors[-1] ], drop = TRUE ))
      x.factor <- lsm[[factors[1]]]
    }
    else { # list of two vectors of factor names
      group <- as.character(interaction( lsm[, factors[[2]] ], drop = TRUE ))
      x.factor <- interaction( lsm[, factors[[1]] ], drop = TRUE )
      xlabs <- paste(factors[[1]], collapse=".")
    }
  }
  x = x.factor
  ## set up horizontal axis
  if( is.xyplot ) {
    x = unfactor( x )
  }
  if( is.factor( x )) {
    xlims = c( 0.5, 0.5 + length( levels( x )))
    if ( sort.mean ) {
      ## factor levels sorted by mean
      o = if (is.null(group))
        order(-lsm$pred)
      else
        order(-tapply(lsm$pred, x.factor, mean))
      lsm = lsm[o,]
      x = ordered( x[o], levels( x )[ o ] )
    }
  }
  else
    xlims = range( pretty( x ))
  width = width * diff( xlims )
  
  ylims = range( pretty( ylims ))

  ## prepare plot panel
  switch( bar.plot, 
         none = {
           mypanel = function(x,y,...) {
             panel.superpose(x,y,...)
             myaxis(x,y,...)
             panelf(x,y,...)
           }
         }, 
         lsd = {
           xposn = xlims[2] - 0.2 * diff( xlims )
           yposn = ( ylim[1] + sum( ylim )) / 3
           mypanel = function(x,y,...) {
             panel.superpose(x,y,...)
             myaxis(x,y,...)
             panelf(x,y,...)
             lsd.bar.lm( object, data, factors, lsm = lsm, 
                        width = width, level = level, 
                        rdf = rdf, xpos=xpos, ypos=ypos, ... )
           }
         }, 
         test =, ci = {
           xx <- unfactor(x, "ordered")
           group <- if (length(factors) == 2) lsm[[factors[2]]] else NULL
           offset = offset * width
           if (!is.null(group) & offset > 0)
             xx <- xx + offset *
               (unclass(group) - (1+length(levels(group)))/2)
           mypanel = function(x,y,...) {
             panel.superpose(x,y,...)
             myaxis(x,y,...)
             panelf(x,y,...)
             se.bar(xx, y, 2*ci, width=width, cap="", ...)
           }
         }, 
         ellipse = {
           ## set up ellipses; adjust ylims for ellipses
           set.ellipse = function( x, object, lsm, rdf, level ) {
             if( is.null( lsm$pdiff )) {
               plet <- pletter( pdiff( object, lsm, rdf = rdf ),
                               level = level )
               plet <- plet[ order( order( - lsm$pred )) ]
             }
             else
               plet <- lsm$pdiff
             
             xx <- unfactor(x, "ordered")
             ux <- unique( xx )
             ellipx <- list()
             o <- order( - lsm$pred )
             for( i in seq( length( ux ))) {
               ii <- xx == ux[i]
               ## there is some problem here with small sets
               pset <- pletter.list( plet[o], subset = order( o )[ii] )
               ellipx[[i]] <- pvalue.ellipse( plet[ii], lsm$pred[ii], 
                                             lsm$se[ii], rdf, 
                                             level = level, pset = pset )
             }
             list( xx = xx, ellipx = ellipx )
           }
           ellipx = set.ellipse( x, object, lsm, rdf, level )
           xx = ellipx$xx
           ellipx = ellipx$ellipx
           for( i in seq( ellipx ))
             if( length( ellipx[[i]] ))
               ylims <- range( ylims, ellipx[[i]][1, ] - ellipx[[i]][2, ],
                              ellipx[[i]][1, ] + ellipx[[i]][2, ] ) 
           ux <- unique( xx )
           width = 3 * width
           tmpfn <- function( y, origin, width ) {
             if( !is.na( y[2] ) & y[2] > 0 )
               panel.ellipse( origin[1], origin[2] + y[1],
                             width, y[2], center = 0, fine = fine )
             return(0)
           }
           mypanel = function(x,y,...) {
             panel.superpose(x,y,...)
             myaxis(x,y,...)
             panelf(x,y,...)
             for( i in seq( ellipx )) {
               tmp <- ellipx[[i]]
               if( length( tmp ) & !is.na( tmp[1] ))
                 apply( tmp, 2, tmpfn, c(ux[i], 0), width )
             }
           }
         })
  if( is.null( group ))
    group = rep( "o", length( x ))
  tmpdata = data.frame( x = x, y = lsm$pred, group = factor( group ))
  ## lattice plot
  if( white )
    trellis.par.set(theme=col.whitebg(), warn = FALSE)
  plottype = if( is.xyplot )
    xyplot
  else
    bwplot
  print( plottype( y ~ x, tmpdata, groups = group,
                  pch = levels( tmpdata$group ),
                  type = type, cex = cex,
                  panel = mypanel,
                  panel.groups = "panel.xyplot",
                  xlab = xlab, ylab = ylab,
                  ylim = ylim, ... ), ... )
  invisible()
}
#####################################################################
margin.plot <- function(object, ...) UseMethod("margin.plot")

margin.plot.default <- function(object, response, group, 
                                orient = "default", average = "full", 
                                xlab=xlabs, ylab=ylabs, ...)
{
  xlabs <- deparse(substitute(object))
  xlabs[2] <- deparse(substitute(group))
  ylabs <- deparse(substitute(response))

  data <- data.frame(response, object, group)
  names(data) <- c(ylabs, xlabs)

  formula.y <- formula( paste(ylabs, ".", sep="~"))
  reduced <- lm( formula.y, data)
  full <- update(reduced, ~.^2, singular.ok=TRUE)
 
  if(is.na(pmatch(average, "full")))
    lsm <- lsmean(reduced, data)
  else
    lsm <- lsmean(full, data)

  orient <- pmatch(orient, c("switch", "both"), nomatch=0)
  if (orient == 1)
    ylabs <- c("", ylabs)
  if(orient != 1)
    ret <- margin.plot.lm(full, reduced, data, factors=xlabs, lsm=lsm, 
                          xlab=xlab[1], ylab=ylab[1], ...)
  if(orient > 0)
    ret <- margin.plot.lm(full, reduced, data, factors=rev(xlabs), lsm=lsm, 
                          xlab=xlab[2], ylab=ylab[2], ...)
  invisible(list(ret=ret, reduced=reduced, full=full))
}
margin.plot.lm <- function(object, reduced=object, 
                           data=eval(object$call$data),
                           factors=all.factors(object, data), 
                           lsm=lsmean(object, data, factors),
                           margin=margin.labels, error="full", 
                           bar.plot="lsd", xlab=xlabs, ylab="",
                           ...)
{
  ## allow for multi-factor horizontal axis
  x.factor <- factors[[1]]
  xlabs <- paste(x.factor, collapse=".")
  factors[[1]] <- xlabs
  x <- interaction( lsm[, xlabs ], drop = TRUE )

  ## use margin values if they correspond to levels of first factor
  ## lazy evaluation: margin=margin.labels
  margin.labels <- precision(lsmean(reduced, data, x.factor))
  lsm$xlab = x
  if (margin[1]) {
    names(margin) <- names(margin.labels)
    lsm[[xlabs]] <- ordered(margin[as.character(x)])
  }
  if(is.na(pmatch(error, "full"))) {
    error <- reduced
    lsm$se <- lsmean(reduced, data, factors)$se
  }
  else
    error <- object
  ## actually want to plot the levels of x, not the values
  ## needs some work
  ## only some plotted?
  myaxis = function(x,y,...) {
    d = !duplicated( lsm$xlab )
    panel.axis( "bottom", x[d], lsm$xlab[d], rot = 0, half = FALSE )
    panel.abline( 0, 1, lty = 3 )
  }
  ret <- int.plot(error, data, factors, lsm=lsm, xlab=xlab, ylab=ylab, 
                  bar.plot=bar.plot,
                  xyplot = myaxis, aspect = 1, ...)

  ## add actual factor 1 levels if using marginal mean values 
 # if (margin[1]) {
 #   if (xaxt != "n") {
 #     mtext(names(margin), 1, 1, at=margin)
 #     right <- par("usr")[1]
 #     mtext(xlabs, 1, 1, at=right)
 #     mtext("margin", 1, 0, at=right)
 #   }
 # }
  invisible(ret)
}
#####################################################################
ci.plot <- function(object, ...) UseMethod("ci.plot")

ci.plot.default <- function(object, response, group,
                            xlab=xlabs, ylab=ylabs, bar.plot="ci", ...)
{
  xlabs <- deparse(substitute(object))
  ylabs <- deparse(substitute(response))
  
  if (!is.factor(object))
    object <- factor(object)
  if( missing(group) ) {
    data <- data.frame(response, object)
    names(data) <- c(ylabs, xlabs)
    formula.y <- formula( paste( ylabs, xlabs, sep="~" ))
    fit <- lm(formula.y, data)
  }
  else {
    glabs <- deparse(substitute(group))
    if (!is.factor(group))
      group <- factor(group)
    data <- data.frame(response, object, group)
    names(data) <- c(ylabs, xlabs, glabs)
    formula.y <- formula( paste( ylabs, ".^2", sep="~" ))
    fit <- lm(formula.y, data)
  }
  invisible(int.plot.lm(fit, data, xlab=xlab, ylab=ylab, bar.plot=bar.plot, ...))
}

ci.plot.lm <- function(object, ..., bar.plot="ci")
  invisible(int.plot.lm(object, ..., bar.plot=bar.plot))
ci.plot.lmer <- ci.plot.lme <-
  function(object, data = eval(object@call$data), ..., bar.plot="ci")
  int.plot.lm(object, data, ..., bar.plot=bar.plot)

ci.width <- function(object, data=eval(object$call$data),
                     factors=all.factors(object, data),
                     lsm=lsmean(object, data, factors, ...),
                     rdf=df.resid(object), level=.05, crit=qt(1-level/2, rdf),
                     ...)
{
  if (level > 0 & level < 1) crit * lsm$se else 0
}
