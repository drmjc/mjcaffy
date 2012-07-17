#' normalize Affymetrix arrays
#' 
#' normalize Affymetrix arrays via dChip, rma, mas5
#' 
#' @param cel.data an \code{AffyBatch} object
#' @param \dots arguments passed to sub-calls
#' 
#' @return An object of class \code{ExpressionSet}
#' 
#' @author Mark Cowley, 27 April 2005
#' @export
#' @rdname normalise-affy
#' @importFrom affy expresso
norm.dchip <- function(cel.data, ...) {
    return( expresso(cel.data, normalize.method="invariantset",
                     bg.correct=FALSE,
                     pmcorrect.method="pmonly", summary.method="liwong") )
}

#' @export
#' @rdname normalise-affy
#' @importFrom affy rma
norm.rma <- function(cel.data, ...) {
    return( rma(cel.data,subset=NULL, verbose=TRUE, destructive = TRUE,
                normalize=TRUE, background=TRUE, bgversion=2, ...) )
}

#' @param scale Value at which all arrays will be scaled to.
#' @export
#' @rdname normalise-affy
#' @importFrom affy mas5
norm.mas5.slow <- function(cel.data, scale=500, ...) {
    return( mas5(cel.data, normalize = TRUE, sc=scale, analysis = "absolute", ...) )
}

#' @export
#' @rdname normalise-affy
#' @importFrom simpleaffy justMAS
norm.mas5 <- function(cel.data, scale=500, ...) {
    return( justMAS(cel.data, tgt=scale, scale=T, ...) )
}


## norm.mas4<- function(cel.data, ...) {
##     require(affy)
##     return( expresso(cel.data, normalize.method="invariantset",
##                      bg.correct=FALSE,
##                      pmcorrect.method="pmonly", summary.method="liwong") )
## }
