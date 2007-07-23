#####################################################################
##
## $Id: ellipse.R,v 1.1 2005/03/24 yandell@stat.wisc.edu Exp $
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
pletter.list <- function( plet, condense = TRUE, subset = seq( plet ) )
{
  plet <- as.character( plet )
  nprod <- length( plet )
  x <- list( )
  i <- 1
  li <- letters[i]
  sprod <- seq( nprod )
  pplet <- paste( plet, collapse = "" )
  anychar <- ".*"
  if( is.null( version$language )) {
    if( pmatch( "win", casefold( version$platform )))
      anychar <- ".*"
  }
  else {
    if( version$language != "R" )
      anychar <- ".*"
  }
    
  while( length( grep( paste( anychar, "[", li, "-z]", anychar, sep = "" ),
                      pplet )) > 0 ) {
    tmp <- grep( paste( anychar, li, anychar, sep = "" ), plet )
    if( length( tmp ))
      x[[li]] <- tmp
    i <- i + 1
    li <- letters[i]
  }
  if( length( x ) == 0 )
    x <- list( sprod )
  ## remove sets that do not concern the subset
  if( length( subset ) < length( plet )) {
    subis <- rep( NA, length( subset ))
    names( subis ) <- as.character( sort( subset ))
    for( i in rev( names( x ))) {
      pm <- cumsum( !is.na( match( sort( subset ), x[[i]] )))
      pm <- pm & is.na( subis )
      if( any( pm ))
        subis[pm]  <- i
      else
        x[[i]] <- NULL
    }
    ## trim sets to contain only subset
    for( i in names( x )) {
      pm <- match( x[[i]], subset )
      pm <- pm[ !is.na( pm ) ]
      if( length( pm ) > 1 )
        x[[i]] <- pm
      else
        x[[i]] <- NULL
    }
  }
  ## condense out singletons and subsets
  if( condense ) {
    ## singletons
    lx = length( x )
    if( lx ) for( i in rev( seq( lx )))
      if( length( x[[i]] ) <= 1 )
        x[[i]] <- NULL
    ## subsets
    if( length( x ) > 1 ) {
      i <- 1
      j <- 2
      while( i < length( x )) {
        tmpij <- any( is.na( match( x[[i]], x[[j]] )))
        tmpji <- any( is.na( match( x[[j]], x[[i]] )))
        while( tmpij & tmpji & j < length( x )) {
          j <- j + 1
          tmpij <- any( is.na( match( x[[i]], x[[j]] )))
          tmpji <- any( is.na( match( x[[j]], x[[i]] )))
        }
        if( !tmpij ) {
          x[[i]] <- NULL
          j <- i + 1
        }
        else {
          i <- i + 1
          if( !tmpji ) {
            x[[j]] <- NULL
          }
        }
      }
    }
  }
  x
}
###########################################################################
pletter <- function( pdiff, level = .05, let = letters )
{
  n <- nrow( pdiff )
  pdiff[ seq( from = 1, to = n*n, by = n+1 )] <- 1
  for( i in seq( n - 1 ))
    pdiff[ i, seq( i + 1, n ) ] <- 0
  
  pdiff <- pdiff > level
  pdiff[ is.na( pdiff ) ] <- TRUE
  prev <- pdiff[,1]
  
  code <- character( n )
  code[ prev ] <- let[1]
  j <- 1
  for( i in seq( 2, n )) if( any( (!prev) & pdiff[,i] )) {
    j <- j + 1
    prev <- pdiff[,i]
    code[prev] <- paste( code[prev], let[j], sep = "" )
  }
  names(code) <- dimnames( pdiff )[[1]]
  code
}
###########################################################################
pletter.ravel <- function( pdiff, level = .05, let = letters )
{
  ## assume have the following:
  ## levels ordered from largest to smallest
  ## col 1 = i
  ## col 2 = j, j > i
  ## col 3 = pij pvalue
  ## rownames = i.j
  pdiff[[3]] <- pdiff[[3]] > level
  pdiff[[3]][ is.na( pdiff[[3]] ) ] <- TRUE
  plevels <- levels( pdiff[[1]] )
  n <- length( plevels )

  prev <- c( TRUE, as.logical( pdiff[[3]][ seq( n - 1 ) ] ))
  
  code <- character( n )
  names(code) <- plevels
    
  code[ prev ] <- let[1]
  k <- 1
  j <- n - 1
  for( i in seq( 2, n - 1 )) {
    newprev <- c( rep( FALSE, i - 1 ), TRUE,
                 as.logical( pdiff[[3]][ j + seq( n - i ) ] ))
    j <- j + n - i
    prev[i-1] <- FALSE
    if( any(( (!prev) & newprev ) | ( prev & (!newprev) ))) {
      k <- k + 1
      prev <- newprev
      code[prev] <- paste( code[prev], let[k], sep = "" )
    }
  }
  code
}
###########################################################################
pdiff.ravel <- function( fit, lsm = lsmean( fit ),
                        pred = outer( lsm$pred[order.pred],
                          lsm$pred[order.pred], "-" ),
                        se = outer( lsm$se[order.pred], lsm$se[order.pred],
                          function( x, y ) sqrt( x^2+y^2 )),
                        rdf = df.resid( fit ))
{
  order.pred <- order( -lsm$pred )
  names.pred <- names( lsm$pred[order.pred] )
  n <- length( order.pred )
  tmprc <- lower.tri( diag( n ))
  indrc <- col( diag( n ))[tmprc]
  ni <- length( indrc )
  pdiff <- data.frame( col = ordered( names.pred[indrc], names.pred ))
  
  difmat <- matrix( 0, ni, n,
                   dimnames = list( NULL, names( lsm$pred )))
  difmat[ seq( ni ) + ni * ( indrc - 1 ) ] <- 1
  indrc <- row( diag( n ))[tmprc]
  difmat[ seq( ni ) + ni * ( indrc - 1 ) ] <- -1
  
  pdiff$row <- ordered( names.pred[indrc], names.pred )

  tmp <- - abs( pred / se )
  2 * pt( tmp, rdf )
  pdiff$pvalue <- 2 * pt( tmp, rdf )[tmprc]
  pdiff
}
###########################################################################
pdiff <- function( fit, lsm = lsmean( fit ),
                  pred = outer( lsm$pred[order.pred], lsm$pred[order.pred], "-" ),
                  se = outer( lsm$se[order.pred], lsm$se[order.pred],
                    function( x, y ) sqrt( x^2+y^2 )),
                  rdf = df.resid( fit ))
{
  order.pred <- order( -lsm$pred )
  tmp <- - abs( pred / se )
  2 * pt( tmp, rdf )
}
###########################################################################
pvalue.ellipse <- function( plet, lsmeans, se, rdf,
                           level = .05, scale = 0.3, full.width = TRUE,
                           pset = pletter.list( plet ))
{
  level <- 1 - level / 2
  ni <- length( lsmeans )
  ns <- length( pset )
  if( ns == 0 )
    return( matrix( NA, 2, 0 ))
  elstuff <- matrix( NA, 2, ns, dimnames = list( c("x","sx"), 1:ns ))
  if( length( rdf ) == 1 )
    rdf <- array( rdf, ni )

  pnames <- names( lsmeans )
  if( is.null( pnames ))
    pnames <- as.character( seq( lsmeans ))

  tmp <- function( r, se, rdf, level, scale, full.width,
                  tmpfn = function( x, y ) sqrt( x^2 + y^2 ))
  {
    if( length( se ) > 1 )
      {
        se <- outer( se, se, tmpfn )
        se <- max( se[ row( se ) > col( se ) ] )
        radius <- scale * se + diff( r ) / 2
      }
    else
      {
        radius <- se <- 0
      }
    rdf <- min( rdf )
    
    if( full.width )
      c( mean( r ), max( radius, se * qt( level, rdf ) / 2 ))
    else
      c( mean( r ), max( radius ))
  }
  for( i in seq( ns )) {
    same <- pset[[i]]
    elstuff[,i] <- tmp( range( lsmeans[same] ), se[same], rdf[same],
                       level, scale, full.width )
  }
  as.matrix( elstuff )
}
###########################################################################
circle <- function(x=0, y=0, radius=1, sx=radius, sy=radius, ...)
  invisible(ellipse(x, y, sx=sx, sy=sy, ...))

ellipse <- function(x=0, y=0, sx=1, sy=1, rho=0, fine=.005, center=.1,
                    show=TRUE, add=TRUE, 
                    xlab="", ylab="", bty="n", xaxt="n", yaxt="n", type="n",
                    lty=1, pty="s", ...)
{
  if(min(sx) < 0 | min(sy) < 0 | max(abs(rho)) > 1)
    stop("sx, sy or rho out of bounds")
  if(fine <= 0 | fine > 1)
    stop("fine out of bounds")
  n <- max(length(x), length(y), length(sx), length(sy), length(rho))
  x <- array(x, n)
  y <- array(y, n)
  lty <- array(lty, n)
  sx <- array(sx, n)
  sy <- array(sy, n)
  rho <- array(rho, n)
  center <- array(center, n)
  prin1 <- sqrt((1+rho)/2)
  prin2 <- sqrt((1-rho)/2)
  xarc <- sin(pi * seq(0, .5, by=4*fine))
  yarc <- c(rev(xarc), -xarc)
  xarc <- c(xarc, yarc, -rev(xarc))
  yarc <- c(yarc, rev(yarc))
  arc <- list()
  for (i in 1:n)
    arc[[i]] <- data.frame(x=x[i]+sx[i]*(prin1[i]*xarc+prin2[i]*yarc), 
                           y=y[i]+sy[i]*(prin1[i]*xarc-prin2[i]*yarc))

  if(show) {
    if(!add) {
      tmp <- options( warn = -1 )
      plot(c(min(x-sx), max(x+sx)), c(min(y-sy), max(y+sy)), 
           bty=bty, xaxt=xaxt, yaxt=yaxt, type=type, xlab=xlab, ylab=ylab,
           pty=pty, ...)
      options( tmp )
    }
    err <- -1
    for (i in 1:n) {
      lines(arc[[i]]$x, arc[[i]]$y, lty=lty[i], err=err, ...)
      if(center[i] > 0) {
        lines(x[i]+center[i]*c(-1, 1)*sx[i]*prin1[i], 
              y[i]+center[i]*c(-1, 1)*sy[i]*prin1[i], err=err)
        lines(x[i]+center[i]*c(-1, 1)*sx[i]*prin2[i], 
              y[i]+center[i]*c(1, -1)*sy[i]*prin2[i], err=err)
      }
    }
  }
  invisible(arc)
}
panel.ellipse <- function(x=0, y=0, sx=1, sy=1, rho=0, fine=.005, center=.1,
                    xlab="", ylab="", bty="n", xaxt="n", yaxt="n", type="n",
                    lty=1, pty="s", ...)
{
  if(min(sx) < 0 | min(sy) < 0 | max(abs(rho)) > 1)
    stop("sx, sy or rho out of bounds")
  if(fine <= 0 | fine > 1)
    stop("fine out of bounds")
  n <- max(length(x), length(y), length(sx), length(sy), length(rho))
  x <- array(x, n)
  y <- array(y, n)
  lty <- array(lty, n)
  sx <- array(sx, n)
  sy <- array(sy, n)
  rho <- array(rho, n)
  center <- array(center, n)
  prin1 <- sqrt((1+rho)/2)
  prin2 <- sqrt((1-rho)/2)
  xarc <- sin(pi * seq(0, .5, by=4*fine))
  yarc <- c(rev(xarc), -xarc)
  xarc <- c(xarc, yarc, -rev(xarc))
  yarc <- c(yarc, rev(yarc))
  arc <- list()
  for (i in 1:n)
    arc[[i]] <- data.frame(x=x[i]+sx[i]*(prin1[i]*xarc+prin2[i]*yarc), 
                           y=y[i]+sy[i]*(prin1[i]*xarc-prin2[i]*yarc))
  err <- -1
  for (i in 1:n) {
    panel.lines(arc[[i]]$x, arc[[i]]$y, lty=lty[i], err=err, ...)
    if(center[i] > 0) {
      panel.lines(x[i]+center[i]*c(-1, 1)*sx[i]*prin1[i], 
                  y[i]+center[i]*c(-1, 1)*sy[i]*prin1[i], err=err)
      panel.lines(x[i]+center[i]*c(-1, 1)*sx[i]*prin2[i], 
                  y[i]+center[i]*c(1, -1)*sy[i]*prin2[i], err=err)
    }
  }
  invisible(arc)
}
norm.ellipse <- function(x, y, level=.05, 
                         radius=sqrt(qchisq(1-level, 2)), ...) {
  mx <- mean(x)
  my <- mean(y)
  sx <- std.dev(x)
  sy <- std.dev(y)
  rho <- if(length(x) > 1) cor(x, y, trim=.00001) else 1
  ellipse(mx, my, sx*radius, sy*radius, rho, ...)
  invisible(list(x=c(mx, sx), y=c(my, sy), rho=rho))
}
