#' plot a distribution of RMA data, split by the affy P/M/A call.
#' 
#' either supply a vector of celfile names as full paths, or
#' supply 2 data.frames (not eSet's) of RMA and MAS5 calls.
#' 
#' @param rma a \code{data.frame} of normalised data. If \code{NULL} then specify \code{celfiles}
#' @param calls a \code{data.frame} of detection calls.  If \code{NULL} then specify \code{celfiles}
#' @param celfiles a vector of cel files.
#' 
#' @return nonw. a plot is created
#' 
#' @author Mark Cowley, 1/4/08
#' @export
#' @importFrom affy read.affybatch justRMA mas5calls
#' @importFrom Biobase exprs
plot_rma_vs_calls <- function(rma=NULL, calls=NULL, celfiles=NULL) {
	
	if( is.null(rma) && !is.null(celfiles) )
		rma <- exprs(justRMA(filenames=basename(celfiles), celfile.path=dirname(celfiles[1])))
	stopifnot( !is.null(rma) )
	
	if( is.null(calls) && !is.null(celfiles)) {
		cel <- read.affybatch(filenames=celfiles)
		calls <- exprs(mas5calls(cel))
	}
	stopifnot( is.null(calls) )
	
	for(i in 1:ncol(rma)) {
		A <- rma[calls[,i]=="A",i]
		M <- rma[calls[,i]=="M",i]
		P <- rma[calls[,i]=="P",i]
		plot( density(rma[,i]), main=colnames(rma)[i], ylim=c(0,0.3), xlab="Expression Level (log2)")
		lines(density(A), col="darkgrey", lty="dashed", lwd=2)
		lines(density(M), col="darkgrey", lty="dotted", lwd=2)
		lines(density(P), col="darkgrey", lty="solid", lwd=2)
		legend("topright", c("all", "A", "M", "P"), col="grey", lwd=c(1,2,2,2),lty=c("solid", "dashed", "dotted", "solid"), bty="none", inset=0.01)
	}

}
