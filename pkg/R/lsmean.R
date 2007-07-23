#####################################################################
##
## $Id: lsmean.R,v 1.1 2005/03/24 yandell@stat.wisc.edu Exp $
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
### fix kludges in lsmean
##################################################################
lsmean <- function(object, ...) UseMethod("lsmean")

lsmean.default <- function(object, ..., factors=all.factors(fit, data)[1:2], 
                           effects=FALSE, se.fit=TRUE, adjust.covar=TRUE)
{
  data <- eval(substitute(data.frame(object, ...)))
  object <- deparse(substitute(object))
  fit <- lm(object~., data, singular.ok=TRUE)
  lsmean.lm(fit, data, factors, effects=effects, se.fit=se.fit,
            adjust.covar=adjust.covar)
}

lsmean.lm <- function( object, data = eval( object$call$data ), 
                      factors = names( predictors[predictors] ), 
                      expr = formula( object ), 
                      contrast = object$contrasts, effects = FALSE, 
                      se.fit = TRUE, adjust.covar = TRUE, 
                      pdiff = FALSE, reorder = FALSE, lsd = FALSE, level = .05, 
                      rdf = df.resid( object ), 
                      coef = coefficients( object ), 
                      cov = vcov( object ), ...)
{
  expr.labels = attr( terms( object ), "term.labels")
  
  ## get predictors from right side of equation
  predictors <- unlist(lapply(data, is.factor))[
    all.names(expr, unique = TRUE, functions = FALSE)[-1]]
  predictors <- predictors[!is.na(charmatch(names(predictors), expr.labels))]

  # set up factors of interest, auxiliary factors and covariates
  covars <- names(predictors[!predictors])
  if( !length( covars ))
    covars <- NULL
  aux.factors <- predictors[predictors]
  if (length(aux.factors) & !is.null( factors )) {
    pm <- pmatch(unlist(factors), names(aux.factors))
    pm <- pm[!is.na(pm)]
  }
  else
    pm <- numeric(0)
  ## R contrasts are sometimes null or just the names
  if( is.null( contrast )) {
    for( i in names( predictors )) {
      if( is.factor( data[[i]] ))
        contrast[[i]] = contrasts( data[[i]] )
    }
  }
  else {
    for( i in names( contrast ))
      if( pmatch( "contr", as.character( contrast[[i]][1] ), nomatch = 0 ))
        contrast[[i]] <- contrasts( data[[i]] )
  }
  if(length(pm)==0) {
    factors <- NULL
    aux.factors <- names(aux.factors)
  }
  else {
    factors <- names(aux.factors[pm])
    aux.factors <- names(aux.factors[-pm])
    
    ## check contrasts if using effects option
    if (effects) {
      for (i in factors) if (any(abs(apply(contrast[[i]], 2, sum)>1e-12)))
        stop(paste("effects=TRUE not implemented for factors", 
                   "without sum-to-zero contrasts (", i, ")"))
    }
  }
  # separate auxiliary factors with respect to sum-to-zero contrasts
  tmp <- length( aux.factors )
  cont.factors <- rep( FALSE, tmp )
  if( tmp ) for (i in seq( tmp ))
    cont.factors[i] <- any(apply(contrasts(data[[aux.factors[i]]]), 2, sum))
  cont.factors <- aux.factors[cont.factors]
  if(length(cont.factors)) if(!is.null(aux.factors))
    aux.factors <- aux.factors[-pmatch(aux.factors, cont.factors, nomatch=0)]

  ## ignore aux factors whose contrasts sum to zero
  fac <- attr( terms( object ), "factors")
  tmp <- c(factors, covars, cont.factors)
  term.factors <- fac[tmp, ]
  if (length(tmp) > 1)
    term.factors <- apply(term.factors, 2, any)
  not.ignore <- 1 - fac[aux.factors, ]
  if (length(aux.factors) > 1)
    not.ignore <- apply(not.ignore, 2, all)
  if (length(aux.factors) & length(not.ignore))
    term.factors <- term.factors & not.ignore
  term.factors <- expr.labels[term.factors]

  ## simplify design: unique factor combinations, covariates at mean level
  ## unique factor combinations
  if(is.null(factors)) {
    newdata <- data[1, ]
    int.factors <- "mean"
    row.names(newdata) <- int.factors
  }
  else {
    int.factors <- as.character(interaction( data[, factors], drop = TRUE ))
    design <- !duplicated(int.factors)
    newdata <- data[design, ]
  }
  # average covariates to mean value
  if (adjust.covar | is.null(factors)) {
    for (i in covars)
      newdata[[i]] <- rep( mean( data[[i]], na.rm = TRUE ), nrow( newdata ))
  }
  else {
    for (i in covars)
      newdata[[i]] <- c( tapply( data[[i]], int.factors, mean ) )
  }
  if(!is.null(factors)) {
    for (i in rev(factors))
      newdata <- newdata[order(newdata[[i]]), ]
    row.names(newdata) <- int.factors <-
      as.character(interaction( newdata[, factors], drop = TRUE ))
  }
  if(!is.null(aux.factors)) {
    for (i in aux.factors )
      levels( newdata[[i]] ) <- levels( data[[i]] )
  }

  ## THIS BREAKS FOR NESTED MODELS
  ## DOES NOT PROPERLY SUM OVER NESTED LEVELS (E.G. IN SPLIT PLOT)
  ## expand matrix for auxiliary factors without sum-to-zero contrasts
  if (!is.null(cont.factors) & length( cont.factors )) {
    factornames <- list()
    for( i in names( data ))
      factornames[[i]] <- levels( data[[i]] )
    cont.levels <- expand.grid( factornames[ cont.factors ])
    cont.n <- nrow(cont.levels)
    newdata <- newdata[rep(int.factors, cont.n), ]
    tmp <- rep(seq(len=cont.n), rep(length(int.factors), cont.n))
    newdata[, cont.factors] <- cont.levels[tmp, ]
  }

  ## get indexes for terms used in LS mean fit
  coef <- coef[ !is.na( coef ) ]
  ## pos <- attr(coef, "assign")
  ## kludge! (R does not have assign as attribute for coef
##  mat <- model.matrix.default( object, newdata, contrast )
  class( expr ) = "formula"
  mat = model.matrix( expr, newdata )
  ## mat <- model.matrix( object )
  pos <- attr( mat, "assign" )
  if( !is.list( pos )) {
    posnames <- attr( terms( object ), "term.labels" )
    posnum <- pos
    pos <- list()
    pos[["(Intercept)"]] <- 1
    
    spos <- seq( length( posnum ))
    for( i in seq( length( posnames )))
      pos[[ posnames[i] ]] <- spos[ i == posnum ]
  }
  index <- if (!effects) pos[["(Intercept)"]]
  for (i in term.factors)
    index <- c(index, pos[[i]])
  # make sure index values are not missing for mat and for coef
  is.coef <- match( dimnames( mat )[[2]], names( coef ), nomatch = 0 )
  mindex <- index[ is.coef[index] > 0 ]
  index <- is.coef[index]

  mat <- t( as.matrix( mat[, mindex] ))
  ## mat <- t(as.matrix(model.matrix.default(expr, newdata, contrast)[, mindex]))
  ## average over levels for auxiliary factors without sum-to-zero contrasts
  if (!is.null(cont.factors) & length( cont.factors )) {
    tmp <- array(mat, c(nrow(mat), ncol(mat)/cont.n, cont.n))
    dimnames(tmp) <- list(dimnames(mat)[[1]], NULL, NULL)
    mat <- apply(tmp, c(1, 2), mean)
  }
  mat <- as.matrix(mat)
  if( nrow( mat ) == 1 )
    mat <- t( mat )
  pred <- coef[index] * mat
  if (length(index) == 1)
    pred <- unlist(pred)
  else
    pred <- apply( pred, 2, sum, na.rm = TRUE )

  ## standard error
  if(se.fit) {
    if(length(index) == 1) 
      se <- sqrt( cov[index, index] ) * abs( mat )
    else
      se <- sqrt( apply( mat * ( cov[index, index] %*% mat ), 2, sum ))
    se.fit <- !any(is.na(se))
  }
  if(is.null(factors) & is.null(covars))
    newdata <- data.frame(pred=matrix(pred, 1, 1), row.names=int.factors)
  else {
    newdata <- as.data.frame(newdata[int.factors, c(factors, covars)])
    dimnames(newdata) <- list(int.factors, c(factors, covars))
    names(pred) <- int.factors
    newdata$pred <- pred
  }
  if(se.fit) {
    if(is.null(factors) & is.null(covars))
      se <- matrix(se, 1, 1)
    newdata$se <- se
    names(newdata$se) <- row.names(newdata)
  }
  n <- length( pred )
  if( pdiff & n > 2 ) {
    tmprc <- lower.tri( diag( n ))
    indrc <- col( diag( n ))[tmprc]
    ni <- length( indrc )
    order.pred <- order( -pred )
    names.pred <- names( pred[order.pred] )
    pdiff <- data.frame( col = ordered( names.pred[indrc], names.pred ))
    difmat <- matrix( 0, ni, n, 
                     dimnames = list( NULL, names( pred )))
    difmat[ seq( ni ) + ni * ( indrc - 1 ) ] <- 1
    indrc <- row( diag( n ))[tmprc]
    difmat[ seq( ni ) + ni * ( indrc - 1 ) ] <- -1
    pdiff$row <- ordered( names.pred[indrc], names.pred )
    
    mat <- mat[, order.pred] %*% t( difmat )
    difpred <- apply( coef[index] * mat, 2, sum, na.rm = TRUE )
    difse <- apply( mat * ( cov[index, index] %*% mat ), 2, sum )
    pdiff$pvalue <- 2 * pt( - abs( difpred / sqrt( difse )), rdf )
    plet <- pletter.ravel( pdiff, level = level )
    newdata$pdiff <- plet[ names( pred ) ]
    if(lsd & !any(is.na(difse)))
      newdata$lsd <- rep( qt( 1 - level / 2, rdf ) * sqrt( mean( difse )), 
                         length( pred ))
  }
  if( reorder )
    newdata <- newdata[ order.pred, ]
 newdata
}
terms.lme = function( x, ... )
{
  ## kludge for lme4 because it lacks terms.lme
  pkg = attr( class( x ), "package" )
  if( is.null( pkg ) | pkg == "nlme" )
    terms( formula( x ))
  else
    slot( x, "terms" )
}
terms.lmer = function( x, ... )
{
  x@terms
}

## nlme or old lme4
lsmean.lme <- function(object, data, factors, ..., 
                       rdf = df.resid(object), 
                       coef = fixef(object),
                       cov = vcov(object))
{
  if( missing( data ) | is.null( data ))
    data = eval( slot(object,"call")$data )

  if( !missing( factors ) & missing( rdf )) {
    rdf <- df.resid( object )
    lsmean.lm( object, data, factors, ..., rdf = rdf, coef = coef, cov = cov )
  }
  else
    lsmean.lm( object, data, ..., rdf = rdf, coef = coef, cov = cov )
}
## new lme4
lsmean.lmer <- function(object, data = eval(object@call$data), factors,
                        expr = terms(object), ..., 
                        rdf = df.resid(object), 
                        coef = fixef(object),
                        cov = as.matrix(vcov(object)))
{
  if( missing( factors ))
    lsmean.lm( object, data,, expr = expr, ...,
              rdf = rdf, coef = coef, cov = cov )
  else
    lsmean.lm( object, data, factors, expr = expr, ...,
              rdf = rdf, coef = coef, cov = cov )
}

lsmean.listof <- function( object, data = eval( attr( object, "call" )$data ), 
                          factors = names( fit.cont ), 
                          stratum = length( object ), 
                          expr = express, 
                          contrast = fit.cont, ... )
{
  express <- attr( object, "call" )$formula
  tmp <- gsub( " . Error(.*)", "", 
              as.character( express[3] ))
  express <- formula( paste( express[2], tmp, sep = "~" ))
  fit.cont <- attr( object, "contrasts" )
  for( i in names( fit.cont ))
    if( pmatch( "contr", as.character( fit.cont[[i]][1] ), nomatch = 0 ))
      fit.cont[[i]] <- contrasts( data[[i]] )
  cat( "WARNING: lsmean does not handle Error() properly. SEs turned off\n" )
  lsmean.lm(object[[stratum]], data=data, factors=factors, expr=expr, 
            contrast=contrast, ..., se.fit = FALSE )
}
