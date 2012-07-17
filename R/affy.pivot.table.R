#' merge expression and calls together
#' 
#' @param exp a \code{matrix}-like object of expression levels
#' @param calls a \code{matrix}-like object of detection calls (P/M/A)
#' @param file the output filename
#' @return none.
#' @author Mark Cowley
#' @export
affy.pivot.table <- function(exp, calls, file=NULL) {
	if( !is.matrix.like(exp) ) exp <- exprs(exp)
	if( !is.matrix.like(calls) ) calls <- exprs(calls)
	exp <- round(exp,4)
	
    res <- interleave.columns(exp, calls)
    colnames(res) <- sub(".CEL", "", colnames(res))
    tmp2 <- rep(c("_Signal", "_Detection"), ncol(exp))
    colnames(res) <- paste(colnames(res), tmp2, sep="")

    if( !is.null(file) ) {
		write.delim(res, file, row.names=TRUE)
	}

    invisible(res)
}
