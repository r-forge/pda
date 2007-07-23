#####################################################################
##
## $Id: source.R,v 1.1 2005/03/24 yandell@stat.wisc.edu Exp $
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
### fix kludges in lsmean
### unfactor, ellipse
##################################################################
arrow <- function(x, y, head=3, tip=.025, width=tip/2, fine=.005, ...)
{
  lines(x, y, ...)
  n <- length(x)
  u <- par()$usr
  diag <- (u[1]-u[2])^2+(u[3]-u[4])^2
  
  if (head%%2)
    dotip(x[1], x[2]-x[1], y[1], y[2]-y[1], diag, tip, width, ...)
  if (head>1)
    dotip(x[n], x[n-1]-x[n], y[n], y[n-1]-y[n], diag, tip, width, ...)
}
dotip <- function(x, dx, y, dy, diag, tip, width, ...)
{
  dd <- sqrt((dx^2+dy^2)/diag)
  a <- tip / dd
  b <- width / (2*dd)
  tx <- x + a*dx + c(1, -1)*b*dy
  ty <- y + a*dy + c(-1, 1)*b*dx
  polygon(c(x, tx), c(y, ty), ...)
}

arc.arrow <- function(diam, head=3, fine=.005, x=diam, y=.5+2*arc,
                      arc=(.5-diam)/4, ...)
{
  tmp <- circle(x, y, diam, show=FALSE, fine=fine)[[1]]
  tmp <- tmp[tmp$y>=y-arc&tmp$y<=y+arc&tmp$x<x, ]
  arrow(tmp$x, tmp$y, head, ...)
  invisible(tmp)
}
path3 <- function(cell=rep("", 4), xlab="", xlim=c(-.1,1), diam=.15, ylab="", ...)
{
  tmp <- diam / sqrt(5)
  
  plot(xlim, 0:1, type="n", axes=FALSE, xlab="", ylab=ylab, ...)
  mtext(xlab, 1)
  
  if(nchar(cell[1])){
    circle(1-diam, .5, diam, cen=0)
    text(1-diam, .5, cell[1])
  }
  if(nchar(cell[2])){
    circle(diam, 1-diam, diam, cen=0)
    text(diam, 1-diam, cell[2])
    arrow(c(diam+2*tmp, 1-diam-2*tmp), c(1-diam-tmp, .5+tmp), 2)
  }
  if(nchar(cell[3])){
    circle(diam, .5, diam, cen=0)
    text(diam, .5, cell[3])
    arrow(c(2*diam, 1-2*diam), c(.5,.5), 2)
  }
  if(nchar(cell[4])){
    circle(diam, diam, diam, cen=0)
    text(diam, diam, cell[4])
    arrow(c(diam+2*tmp, 1-diam-2*tmp), c(diam+tmp, .5-tmp), 2)
  }
}
##bekk <- read.table("../src/bekk.dat", header=TRUE)
##bekk$lab <- factor(bekk$lab)
##bekk$mat <- factor(bekk$mat)
## see help(trans.plot)
box.fig <- function(fit, xpos=1, width=.2, jitter=width)
{
  fit <- boxplot(fit, plot=FALSE)
  if( is.null( fit$stats ))
    fit <- fit[[1]]
  tmp <- fit$stats
  polygon(xpos+width*c(-1, 1, 1, -1), rep(tmp[c(2, 4)], rep(2, 2)))
  lines(xpos+width*c(-1,1), rep(tmp[3],2))
  lines(xpos+width*c(-1,1), rep(tmp[1],2), lty=2)
  lines(xpos+width*c(-1,1), rep(tmp[5],2), lty=2)
  lines(rep(xpos,2), tmp[1:2], lty=2)
  lines(rep(xpos,2), tmp[4:5], lty=2)
  if (!is.null(fit$out) & length( fit$out ))
    points(xpos+runif(length(fit$out), -jitter, jitter), fit$out)
  fit
}
effect.plot <- function(object, effect=effect.fit(object, adjust=adjust, ...),
                        adjust=TRUE, boxplot.min=10, jitter=.1, width=jitter,
                        xlim=xlims, ylim=ylims, xaxt="i",
                        xlab="Terms", ylab=ylabs, ...)
{
  nterms <- length(effect)
  ## eliminate degenerate effects
  for (i in 1:nterms)
    if (all(is.na(effect[[i]])))
      effect[[i]] <- NULL
  nterms <- length(effect)
  
  if (adjust)
    ylabs <- "MS Adjusted Effects"
  else
    ylabs <- "Effects"
  
  xlims <- c(1-max(jitter, width), nterms+max(jitter, width))
  ylims <- range(unlist(effect))
  plot(xlim, ylim, type="n", xaxt="n", xlab=xlab, ylab=ylab, ...)
  if( xaxt != "n" )
    axis(1, 1:nterms, names(effect))
  abline(h=0, lty=2)
  par(lty=1)
  
  for (i in 1:nterms) {
    neff <- length(effect[[i]])
    if (neff < boxplot.min) {
      name1 <- names(effect[[i]])
      if(is.null(name1))
        points(i+runif(neff, -jitter, jitter), effect[[i]])
      else
        text(i+runif(neff, -jitter, jitter), effect[[i]], labels=name1)
      
    }
    else
      box.fig(effect[[i]], i, width, jitter)
  }
  invisible(effect)
}


## effects contrasts

effect.fit <- function(object, ...) UseMethod("effect.fit")

effect.fit.default <- function(object, data=eval(object$call$data),
                               factors=all.factors(object, data), adjust=TRUE,
                               influence=FALSE,
                               terms.fit =
                               as.data.frame(predict(object, type="terms")),
                               ...)
{
  if (adjust)
    tmp <- function(x, f=NULL) {
      if (is.null(f)) {
        nx <- tapply(x, x, length)
        ux <- as.numeric(names(nx))
        ff <- NULL
      }
      else {
        nx <- tapply(x, f, length)
        ux <- tapply(x, f, mean)
        ff <- tapply(f, f, unique)
        if (is.factor(f))
          ff <- levels(f)[ff]
      }
      lx <- length(ux)
      x <- ux * sqrt(nx*lx/(lx-1))
      list(x=x, f=ff)
    }
  else
    tmp <- function(x, f=NULL) {
      if(is.null(f)) {
        x <- unique(x)
        ff <- NULL
      }
      else {
        x <- tapply(x, f, mean)
        ff <- tapply(f, f, unique)
        if (is.factor(f))
          ff <- levels(f)[ff]
      }
      list(x=x, f=ff)
    }
  
  se <- .01*diff(range(terms.fit))
  eff <- list()
  for (i in names(terms.fit)) {
    tmpterm <- tmp(precision(terms.fit[[i]], se), data[[i]])
    eff[[i]] <- tmpterm$x
    if (!is.na(match(i, factors)))
      names(eff[[i]]) <- as.character(tmpterm$f)
    else
      {
        ## for some reason this does not work in R
        names(eff[[i]]) <- NULL
        names( eff[[i]] ) <- rep( "o", length( eff[[i]] ))
      }
  }

  tmp <- resid(object)
  if (influence)
    eff[["resid"]] <- tmp / sqrt(1-lm.influence(object)$hat)
  else {
    if (adjust)
      eff[["resid"]] <- tmp * sqrt(length(tmp)/df.resid(object))
    else
      eff[["resid"]] <- tmp
  }
  eff
}

effect.fit.listof <- function(object, data=eval(attr(object, "call")$data), 
                              factors=all.factors(object, data),
                              terms.fit=fit.proj, stratum, ...)
{
  fit.proj <- as.data.frame(proj(object)[[stratum]])
  fit.proj$Residuals <- NULL
  effect.fit.default(object[[stratum]], data=data, factors=factors, 
                     terms.fit=terms.fit, ...)
}
df.resid <- function(object, ...) UseMethod("df.resid")

df.resid.default <- function(object, ...)
{
  x <- list(object, ...)
  if(length(x) == 1 && is.list(x[[1]]))
    x <- x[[1]]
  x <- tapply(x[[1]], x, length) - 1
  x[is.na(x)] <- 0
  sum(x)
}

df.resid.lm <- function(object, resid=residuals(object), ...)
{
  rdf <- object$df.resid
  if(is.null(rdf)) {
    p <- object$rank
    if(is.null(p))
      p <- sum(!is.na(coefficients(object)))
    rdf <- length(resid) - p
  }
  rdf
}
## nlme and old lme4
df.resid.lme <- function(object, resid = residuals( object ), 
                         factors = attr( terms( formula( object )), 
                           "term.labels" ), ...)
{
  rdf = object$fixDF
  if( !is.null( rdf ))
    rdf = max( terms( rdf ))
  if(is.null(rdf)) {
    p = object$coef
    p = sum( !is.na(
      if( is.null( p$fixed ))
        unlist( coef( object ))
      else
        p$fixed
      ))
    rdf <- length(resid) - p
  }
  rdf
}
## new lme4
df.resid.lmer <- function(object, ...)
{
  ## this is a guess
#  1 + rev(diff(object@nc))[1]
  length(object@y) - attr(logLik(object), "df") - 1
}
std.dev <- function(object)
{
  ## adapted from summary.lm
  sd <- summary(object)$sigma
  if(is.null(sd)) {
    resid <- residuals(object)
    if(is.null(resid))
      sd <- sqrt(var(object))
    else {
      rdf <- df.resid(object)
      if(rdf > 0)
        sd <- sqrt(sum(resid^2)/rdf)
      else
        sd <- NA
    }
  }
  if(!is.na(sd)) sd
  else 0
}
sample.size <- function(object, ...) UseMethod("sample.size")

sample.size.default <- function(object, ...)
{
  x <- list(object, ...)
  if(length(x) == 1 && is.list(x[[1]]))
    x <- x[[1]]
  x <- tapply(x[[1]], x, length)
  x[is.na(x)] <- 0
  x
}

sample.size.lm <- function(object, data=eval(object$call$data), 
                           factors=all.factors(object, data),
                           na.action=na.omit, ...)
{
  n <- length(factors)
  if (n < 1)
    length(data[[model.response(object)]])
  else {
    replications(formula(paste("~", paste(factors, collapse=":"), sep="")), 
                 data, na.action = na.action)[[1]]
  }
}
expt <- function(cell, effects, inter, grand=0, design=c("rcbd", "crd", "bibd"), 
                 reps=4, factors=c("A", "B"), tukey=FALSE, sd=c(1, .5),
                 anova=TRUE, bibd)
{
  if(missing(cell)){
    if(missing(effects))
      stop("must specify cell means or effects")
    if(length(effects)!=2)
      stop("only two factor effects allowed")
    marg <- expand.grid(effects)
    cell <- grand + apply(marg, 1, sum)
    if(!missing(inter))
      cell <- cell + array(c(inter), length(cell))
    else {
      if(tukey!= 0)
        cell <- cell + tukey*apply(marg, 1, prod)
    }
    cell <- array(cell, unlist(lapply(effects, length)))
    tmp <- function(x) {
      n <- names(x)
      if(is.null(n))
        letters[seq(x)]
      else
        n
    }
    dimnames(cell) <- lapply(effects, tmp)
    factors <- names(effects)
  }
  nt <- length(cell)
  if (!missing(bibd)) {
    if (bibd >= nt | bibd < 2)
      stop(paste("bibd =", bibd, "must be between 2 and", nt))
    design <- "bibd"
  }
  
  sd <- array(c(sd, 0), 2)
  des <- pmatch(design[1], c("rcbd", "crd", "bibd"))
  if (des == 1 & reps > 1) {
    ## rcbd
    y <- rep(rnorm(reps, 0, sd[1]), rep(nt, reps)) + rep(cell, reps) +
      rnorm(reps*nt, 0, sd[2])
    block <- "reps +"
  }
  else if (des == 3) {
    perms <- function(x, k) {
      n <- length(x)
      if (n == 1) n <- x
      if (n < k)
        stop("length of x must be greater than k")
      a <- matrix(c(rep(FALSE, n-k), rep(TRUE,k)), n, 1)
      if (k < n & k > 0) for (i in (n-k-1):0) {
        b <- c(rep(FALSE,i), TRUE)
        d <- perms(n-i-1, k-1)
        a <- cbind(a, rbind(matrix(b, i+1, ncol(d)), d))
      }
      if (length(x) == 1)
        a
      else
        matrix(x, n, ncol(a))[a]
    }
    ## bibd
    numib <- choose(nt, bibd)
    y <- rep(perms(cell, bibd), reps)
    reps <- reps * numib
    nt <- bibd
    y <- y + rep(rnorm(reps, 0, sd[1]), rep(nt, reps)) +
      rnorm(reps*nt, 0, sd[2])
    block <- "reps +"
  }
  else {
    ## crd
    y <- rep(cell, reps) + rnorm(reps*nt, 0, sqrt(sd[1]^2+sd[2]^2))
    block <- ""
  }
  data <- data.frame(reps=factor(rep(1:reps, rep(nt, reps))))
  tmp <- list(ncol(cell), rep(nrow(cell), ncol(cell)))
  for (i in 1:2) {
    cellnames <- dimnames(cell)[[i]]
    if (is.null(cellnames))
      cellnames <- as.character(seq(dim(cell)[i]))
    data[[factors[i]]] <- factor(rep(rep(cellnames, tmp[[i]]), reps))
  }
  data$y <- as.numeric(y)
  if(anova) {
    red <- as.formula(paste("y ~", block, factors[1], "+", factors[2]))
    full <- as.formula(paste("y ~", block, factors[1], "*", factors[2]))
    fullfit <- aov(full, data)
    print(summary(fullfit))
    list(full=fullfit, reduced=aov(red, data), data=data)
  }
  else
    data
}
all.factors <- function (object, data=eval(object$call$data))
{
  predictors <- all.names(formula(object), unique=TRUE, functions=FALSE)[-1]
  factors <- unlist(lapply(data, is.factor))[predictors]
  names(factors[factors])
}

all.covars <- all.covariates <- function (object, data=eval(object$call$data))
{
  predictors <- all.names(formula(object), unique=TRUE, functions=FALSE)[-1]
  factors <- unlist(lapply(data, is.factor))[predictors]
  names(factors[!factors])
}

all.predictors <- function (object, data=eval(object$call$data))
{
  all.names(formula(object), unique=TRUE, functions=FALSE)[-1]
}
##model.response <- function (object, data=eval(object$call$data))
##{
##  all.names(formula(object), functions=FALSE)[1]
##}
factor.labels <- function(object, data) {
  p <- pmatch(names(data), labels(object))
  labels(object)[p[!is.na(p)]]
}
unfactor <- function(x.factor, class="factor")
{
  ## adapted from interaction.plot
  ok <- inherits(x.factor, substitute(class))
  if (ok) {
    xval <- as.character( x.factor )
    ## try to interpret as numeric (turning off warning)
    tmp <- options( warn = -1 )
    nval <- as.numeric( xval )
    options( tmp )
    ok <- !all( is.na( nval ))
    if( ok )
      xval <- nval
  }
  if (ok)
    xval
  else
    unclass(x.factor)
}
lsd.bar <- function(object, ...) UseMethod("lsd.bar")

lsd.bar.default <- function(object, response, group, ...)
{
  xlabs <- deparse(substitute(object))
  ylabs <- deparse(substitute(response))

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
  invisible(lsd.bar.lm(fit, data, ...))
}
lsd.bar.lme <- lsd.bar.lmer <- function(object, data = eval(object@call$data), ...)
  lsd.bar.lm(object, data, ...)

lsd.bar.lm <- function(object, data=eval(object$call$data), 
                       factors=all.factors(object, data),
                       lsm=lsmean(object, data, factors, pdiff=TRUE, lsd=TRUE), 
                       xpos=xposn, ypos=yposn, rdf=df.resid(object), 
                       level=.05,
                       crit=if(rdf > 0)
                         qt(1-level/2, rdf)
                       else
                         1, 
                       cap=paste(100*level, "%\nLSD", sep=""), mod=mods,
                       tol=1e-5, ...)
{
  if( crit == 1 ) {
    if( missing( cap )) cap = "SE"
  }
  if( cap == "SD" ) {
    mod = ""
    rdf = 1
    crit = 1
  }
  usr <- par("usr")
  xlen <- usr[2] - usr[1]
  xposn <- usr[2] - .2 * xlen
  yposn <- (2 * usr[3] + usr[4]) / 3

  if( is.null( lsm$se ))
    lsm$se = rep( 0, nrow( lsm ))
  if( all( is.null( lsm$se )))
    lsm$se = rep( 0, nrow( lsm ))
  
  if( is.null( lsm$lsd )) {
    lsd <- sqrt(mean(lsm$se[lsm$se>0]^2))
    lsd <- crit * sqrt(2) * lsd
  }
  else {
    lsd <- lsm$lsd[1]
    mods <- "*"
  }
  if( !is.null( lsm$se ))
    mods <- if( var( lsm$se, na.rm = TRUE ) > tol * max( abs( lsm$pred )))
      "*"
    else
      ""

  se.bar(xpos, ypos, lsd, mod=mod, cap=cap, ...)
  invisible(list(lsd=lsd, rdf=rdf, level=level))
}
nested <- function(object, ...) UseMethod("nested")

nested.default <- function( object, data = eval( attr( object, "call" )$data ), 
                           project = proj( object ), ... )
  nested.aovprojlist( project, data, ... )

## this assumes that the right list of factors goes along with the nesting

nested.aovprojlist <- function( object, data, factors, response="y", 
                               nesting=c(1, ncol(ref)), ref = projref( object ),
                               ... )
{
  newdata <- data[, factors]
  newdata[[response]] <- apply( ref[ , nesting ], 1, sum )
  projwt( newdata, factors, ... )
}

projwt <- function( data, factors,
                   u = interaction( data[, factors], drop = TRUE ), 
                   weight = 1, covars )
{
  weight <- array( weight, length( u ) )
  newdata <- as.data.frame( data[ !duplicated( u ), ] )
  if ( !missing( covars )) for (i in covars)
    newdata[[i]] <- c(tapply( data[[i]], u, mean ))
  newdata$weight <- c( tapply( weight, u, sum ))
  row.names( newdata ) <- unique( u )
  newdata
}

projref <- function( project )
{
  ref <- data.frame(Reference=c(project[[1]]))
  for (i in names(project)[-1])
    ref[[i]] <- apply( project[[i]], 1, sum )
  ref
}
nested.tree <- function( object, nesting=names(means), 
                        means=nested.means(object, nesting[-1]), offset = 0.05,
                        small = 0.02, 
                        estimate="estimate", bar="stderr", jitter=offset/2, ...)
{
  nn <- length(nesting)

  plot( c(0, 1+nn), range( means[[nn]][[estimate]] ), 
       type = "n", xaxt = "n", ... )
  mtext(nesting, 1, 0, at = seq(nn) )

  small <- rep(small, nn)
  offset <- rep(offset, nn)
  jitter <- rep(jitter, nn)

  mean2 <- means[[1]] 
  estim2 <- mean2[[estimate]]
  for ( i in seq(nn-1) ) {
    mean1 <- mean2
    mean2 <- means[[i+1]]
    estim1 <- estim2
    estim2 <- mean2[[estimate]]
    nester <- nesting[i]
    nest2 <- mean2[[nester]]
    nest1 <- mean1[[nester]]
    nr <- nrow(mean1)
    x <- rep(i, nr)
    text( x-offset[i], uncollide(estim1, small[i]), as.character(nest1), adj = 1 )
    for ( j in seq(nr) )
      for (k in estim2[ nest2==nest1[j] ])
        lines(seq(i, length=2), c(estim1[j], k))
    points( x, estim1 )
    if (jitter[i] > 0)
      x <- x + runif( nr, -jitter[i], jitter[i] )
    se.bar( x, estim1, 2 * mean1[[bar]], cap = "" )
  }
  nester <- nesting[nn]
  invisible(text( nn + offset[nn], uncollide( estim2, small[nn] ), 
                 as.character( mean2[[nester]] ), adj = 0 ))
}

nested.means <- function(object, nesting, parameter="parameter", 
                         carry=c("estimate", "stderr", "ddf"), ...)
{
  x <- list( object[ 1, carry ] )
  root <- as.character(object[ 1, parameter ])
  x[[1]][[root]] <- names(x) <- root
  for (i in seq(nesting)) {
    x[[nesting[i]]] <- object[ object[[parameter]] == nesting[i], 
                              c(nesting[1:i], carry) ]
  }
  x[[2]][[root]] <- rep(root, nrow(x[[2]]))
  x
}

power.curve <- function(n = 100, sigma = 1, alpha = 0.05,
                        dstd = seq(-5, 5, length = 100), df)
{
  if(missing(df)) {
    cstd <- qnorm(1 - alpha/2)
    data.frame(dif = (dstd * sigma)/sqrt(n/2), 
               pow = 1 - pnorm(cstd - dstd) + pnorm( - cstd - dstd))
  }
  else {
    cstd <- qt(1 - alpha/2, df)
    data.frame(dif = (dstd * sigma)/sqrt(n/2), pow = 1 - pt(cstd - dstd, df) +
               pt( - cstd - dstd, df))
  }
}
precision <- function(x, se, digits = 0)
{
  if(is.list(x)) {
    se <- x$se
    p <- x$pred
    if(is.data.frame(x))
      names(p) <- names(se) <- dimnames(x)[[1]]
    x <- p
  }
  if(length(se) == 1 && se == 0) {
    x[abs(x) < 1e-05] <- 0
  }
  else {
    if(length(se) > 1 && length(se) != length(x))
      stop("x and se must have the same length")
    scale <- 10^(digits + 1 - round(log10(se)))
    p <- round(scale * x)/scale
    p[is.na(p)] <- x[is.na(p)]
    p[abs(p) < 1e-05] <- 0
    names(p) <- names(x)
    p
  }
}

##precision <- function(value, se)
##{
##  shifts <- 10^(1-round(log10(2*se)))
##  round(value*shifts) / shifts
##}
  
se.bar <- function(xpos, ypos, bar, mod="", cap="SE", adj=0, width=.01, lty=1,
                   horiz=FALSE, ...)
{
  if(!all(is.na(bar)) && max(bar) > 0 && min(bar) >= 0) {
    if(length(xpos) != length(ypos))
      stop("xpos and ypos must have the same length")
    if(length(xpos) != 1 && length(bar) != 1 && length(bar) > length(xpos))
      stop("xpos and bar must have the same length or one must be scalar")

    bar <- array(bar, length(xpos))
    xpos <- xpos[!is.na(bar)]
    ypos <- ypos[!is.na(bar)]
    bar <- bar[!is.na(bar)]
    lty <- rep( lty, length( bar ) )
    if (horiz) {
      top <- xpos+bar*.5
      bot <- xpos-bar*.5
      for (i in 1:length(ypos))
        panel.lines(c(bot[i], top[i]), rep(ypos[i], 2))
      
      left <- right <- ypos
      if (width > 0) {
        left <- ypos-width
        right <- ypos+width
        for (i in 1:length(ypos)) {
          panel.lines(rep(bot[i], 2), c(left[i], right[i]))
          panel.lines(rep(top[i], 2), c(left[i], right[i]))
        }
      }
    }
    else {
      top <- ypos+bar*.5
      bot <- ypos-bar*.5
      for (i in 1:length(xpos))
        panel.lines(rep(xpos[i], 2), c(bot[i], top[i]), lty=lty[i])
      
      left <- right <- xpos
      if (width > 0) {
        width <- width
        left <- xpos-width
        right <- xpos+width
        for (i in 1:length(xpos)) {
          panel.lines(c(left[i], right[i]), rep(bot[i], 2))
          panel.lines(c(left[i], right[i]), rep(top[i], 2))
        }
      }
      if (cap != "") {
        if(adj)
          panel.text(left, ypos, paste(cap, mod, "=", sep=""), adj=1)
        else
          panel.text(right, ypos, paste("=", cap, mod, sep=""), adj=0)
      }
    }
  }
  list(xpos=xpos, ypos=ypos, bar=bar, mod=mod, cap=cap, adj=adj, width=width)
}
burst <- function(root, branch, pos=c(2, 1))
{
  for (i in branch)
    lines(c(root, i), pos)
}
branch <- function(root, limb, pos=1:2, flip=FALSE)
{
  if (flip)
    for (i in limb)
      lines(pos, c(root, i))
  else
    for (i in limb)
      lines(c(root, i), pos)
}
splitplot <- function(response, nest.factors, covars, weights, data)
{
  nest.int <- interaction( data[, nest.factors ], drop = TRUE )
  ## need to check that nest.factors are factors
  ## how to add weighting to fit?
  formula.y <- formula( paste( response, "~", 
                              paste( nest.factors, sep = "*" )))
  nest.fit <- lm( formula.y, data)
  ## whole plot
  Within <- data
  if (missing(weights))
    weights <- "weights"
  if(is.na(pmatch(weights, names(data))))
    Within[[weights]] <- rep(1, nrow(data))
  Within[[response]] <- resid(nest.fit)
  ## split plot
  Between <- data.frame(data[!duplicated(nest.int), ], 
                        row.names=unique(nest.int))
  Between[[response]] <- predict(nest.fit)[!duplicated(nest.int)]

  Between[[weights]] <- tapply(Within[[weights]], nest.int, sum)

  ## covariates
  if(!missing(covars)) {
    covar.data <- data
    for (i in covars) {
      covar.data[[response]] <- data[[i]]
      Within[[i]] <- resid(nest.fit, covar.data)
      Between[[i]] <- predict(nest.fit, covar.data)
    }
  }
  list(Between=Between, Within=Within)
}
subsamp <- function(variance=c(1, 4), diff=1, alpha=.05, 
                    sub=10^seq(0, log10(50), length=100), eu=tot/sub, tot=100)
{
  ems <- variance[1]+variance[2]*sub
  delta <- eu/(2*ems)
  crit <- qt(1-alpha/2, eu-1)
  pow <- pt(abs(diff)*sqrt(delta)+crit, eu-1)-pt(abs(diff)*sqrt(delta)-crit, eu-1)
  crit <- qt(1-alpha/2, eu-1)/sqrt(delta)
  delta <- delta*diff^2
  data.frame(sub=sub, eu=eu, df=eu-1, ems=ems, delta=delta, pow=1-pow)
}

power.curve <- function(n=100, sigma=1, alpha=.05, dstd=seq(-5, 5, length=100))
{
  cstd <- qnorm(1-alpha/2)
  data.frame(dif=dstd*sigma/sqrt(n/2), pow=1-pnorm(cstd-dstd)+pnorm(-cstd-dstd))
}
trans.plot <- function(fit, 
                       predictors=all.predictors(fit)[1:2], 
                       terms=predict(fit, type="terms"), 
                       xlab="product of margins", ylab="interaction", ...)
{
  a <- terms[, predictors[1]]
  b <- terms[, predictors[2]]
  ab <- terms[, paste(predictors[1:2], collapse=":")]
  if (is.null(ab))
    stop(paste(paste(predictors[1:2], collapse=":"), 
               "is not in list of terms"))
  m <- attributes(terms)$constant
  add <- a*b/m
  plot(add, ab, xlab=xlab, ylab=ylab, ...)
  data <- data.matrix(list(x=add, y=ab))
  fit <- lm(y~x, data)
  lines(sort(add), predict(fit)[order(add)])
  c(power=1-coef(fit)["x"], se=sqrt(vcov(fit)["x","x"]))
}
tukey.plot <- function(x.factor, response, group, ..., orient="default", 
                       xlab=xlabs, ylab=ylabs)
{
  xlabs <- deparse(substitute(x.factor))
  xlabs[2] <- deparse(substitute(group))
  ylabs <- c(deparse(substitute(response)), "")

  data <- data.frame(y=response, a=x.factor, b=group)
  names(data) <- c(ylabs[1], xlabs)

  reduced <- lm(y~., data)
  data$inter <- predict(reduced)^2
  formula.y <- formula( paste( paste( ylabs[1], xlabs[1], sep = "~" ), 
                              xlabs[2], "inter", sep = "+" ))
  full <- lm( formula.y, data, singular.ok=TRUE)
  lsm <- lsmean(full, data, xlabs, adjust.covar=FALSE)

  orient <- pmatch(orient, c("switch", "both"), nomatch=0)
  if(orient != 1)
    ret <- margin.plot.lm(full, reduced, data, factors=xlabs, lsm=lsm, 
                          xlab=xlab[1], ylab=ylab[1], ...)
  if(orient > 0)
    ret <- margin.plot.lm(full, reduced, data, factors=rev(xlabs), lsm=lsm, 
                          xlab=xlab[2], ylab=ylab[2], ...)
  invisible(list(reduced=reduced, full=full, data=data, lsm=lsm, ret=ret))
}

mandel.plot <- function(x.factor, response, group, ..., orient="default", 
                        xlab=xlabs, ylab=ylabs)
{
  xlabs <- deparse(substitute(x.factor))
  xlabs[2] <- deparse(substitute(group))
  ylabs <- c(deparse(substitute(response)), "")

  data <- data.frame(y=response, a=x.factor, b=group)
  names(data) <- c(ylabs[1], xlabs)

  reduced <- lm(y~., data)
  data$inter <- unfactor(lsmean(reduced, data,
                                se.fit=FALSE)[[1]][codes(x.factor)])
  formula.y <- formula( paste( paste( ylabs[1], xlabs[1], sep = "~" ), 
                              xlabs[2], paste( xlabs[2], "inter", sep = ":"), 
                              sep = "+" ))
  full <- lm( formula.y, data, singular.ok = TRUE )
  lsm <- lsmean(full, data, xlabs, adjust.covar=FALSE)

  orient <- pmatch(orient, c("switch", "both"), nomatch=0)
  if(orient != 1)
    ret <- margin.plot.lm(full, reduced, data, factors=xlabs, lsm=lsm, 
                          xlab=xlab[1], ylab=ylab[1], ...)
  if(orient > 0)
    ret <- margin.plot.lm(full, reduced, data, factors=rev(xlabs), lsm=lsm, 
                          xlab=xlab[2], ylab=ylab[2], ...)
  invisible(list(reduced=reduced, full=full, data=data, lsm=lsm, ret=ret))
}
uncollide <- function(x, small = 0.02)
{
  if(length(x) > 1) {
    o <- order(x)
    n <- names(x)
    x <- sort(x)
    d <- (small * diff(range(x)) - diff(x))/2
    d <- d * (d > 0)
    x[o] <- (x - c(d, 0) + c(0, d))
    names(x) <- n
  }
  x
}
