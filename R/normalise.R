## wrapper to affy's dchip normalisation function
##
## Mark Cowley, 27 April 2005
##
norm.dchip <- function(cel.data, ...) {
    require(affy)
    return( expresso(cel.data, normalize.method="invariantset",
                     bg.correct=FALSE,
                     pmcorrect.method="pmonly", summary.method="liwong") )
}

## wrapper to affy's RMA normalisation function
##
## Mark Cowley, 27 April 2005
##
norm.rma <- function(cel.data, ...) {
    require(affy)
    return( rma(cel.data,subset=NULL, verbose=TRUE, destructive = TRUE,
                normalize=TRUE, background=TRUE, bgversion=2, ...) )
}

## wrapper to affy's MAS5 normalisation function
##
## Mark Cowley, 27 April 2005
##
norm.mas5.slow <- function(cel.data, scale=500, ...) {
    require(affy)
    return( mas5(cel.data, normalize = TRUE, sc=scale, analysis = "absolute", ...) )
}


## wrapper to simpleaffy's MAS5 normalisation function
## This is coded in C for speed.
##
## Mark Cowley, 27 April 2005
##
norm.mas5 <- function(cel.data, scale=500, ...) {
    require( simpleaffy )
    return( justMAS(cel.data, tgt=scale, scale=T, ...) )
}


## norm.mas4<- function(cel.data, ...) {
##     require(affy)
##     return( expresso(cel.data, normalize.method="invariantset",
##                      bg.correct=FALSE,
##                      pmcorrect.method="pmonly", summary.method="liwong") )
## }
