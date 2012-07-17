#' Plot a histogram of affy P/M/A calls.
#' 
#' more correctly, this generates a stacked bar plot, where the counts of
#' P, M and A are above each other, for each array
#' 
#' @param calls a data.frame, or expression set from mas5calls
#' @param col the three colours corresponding to P, M, and A
#' @param grid logical
#' @param plot logical
#' @param legend logical
#' @param names.arg character vector
#' @param las see ?par
#' 
#' @return none; generates a plot
#' 
#' @author Mark Cowley, 10 April 08
#' 
#' @examples
#' \dontrun{
#'   par(mar=c(12,5,4,2), cex.lab=0.8)
#'   hist_affycalls(calls)
#'   grid()
#'   dev.off()
#' }
#' @export
#' @importFrom affy exprs
#' @importFrom mjcgraphics axis.percentiles hgrid
hist_affycalls <- function(calls, col=c("green", "orange", "red"), 
							grid=TRUE, plot=TRUE, legend=TRUE, names.arg=colnames(calls), las=2) {
	
	if( class(calls) == "ExpressionSet" ) {
		calls <- exprs(calls)
	}
	
	counts <- apply(calls, 2, table)[c("P", "M", "A"), ]
	
	if( plot ) {
		# opar <- par(no.readonly=TRUE)
		las <- par()$las
		
		prop <- counts / nrow(calls) * 100
		par(las=las)
		barplot(counts, col=col, main="Distribution of Affy Calls", ylab="Number of Present/Marginal/Absent calls", names.arg=names.arg)
		axis.percentiles(max=nrow(calls))
		
		if( grid ) hgrid()
		if( legend ) legend("topright", c("P","M","A"), fill=col, bg="white", inset=0.02)
		
		par(las=las)
		# par(opar)
	}	

	invisible(counts)
}

#' Plot a cumulative histogram of affy P/M/A calls.
#' 
#' undocumented
#' 
#' @param Pcount undocumented
#' @param calls undocumented
#' @param grid logical
#' @param plot logical
#' 
#' @return none; generates a plot
#' 
#' @author Mark Cowley, 10 April 08
#' 
#' @examples
#' \dontrun{
#'   hist_Pcount_cumulative(calls=calls)
#' }
#' @export
#' @importFrom affy exprs
#' @importFrom mjcgraphics axis.percentiles hgrid
hist_Pcount_cumulative <- function(Pcount=NULL, calls=NULL, grid=TRUE, plot=TRUE) {

	if( is.null(calls) && is.null(Pcount) )
		stop("must supply either the mas5calls, or the Pcount.")
	
	if( is.null(Pcount) ) {
		if( class(calls) == "ExpressionSet" ) {
			calls <- exprs(calls)
		}
	
		Pcount <- apply(calls=="P", 1, sum)
	}
	
	PcountCum <- cumsum( rev( table(Pcount) ) ) 
	names(PcountCum)[2:length(PcountCum)] <- paste(">=", names(PcountCum)[2:length(PcountCum)])
	
	if( plot ) {
		opar <- par(no.readonly=TRUE)
		par(las=2)
		barplot(PcountCum, main="Cumulative Pcount", xlab="Number of arrays", ylab="Frequency", col="white")
		axis.percentiles(side=4, length(Pcount))
		if( grid ) hgrid()
		par(opar)
	}

	invisible( PcountCum )
}

