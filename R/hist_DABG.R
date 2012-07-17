#' Plot a histogram of DABG values
#' 
#' Plot a histogram of DABG values, thresholded at 10^-5, 10^-4, 10^-3, 10^-2,
#' 0.1 and 1.0
#' 
#' @param dabg the P-valuse from computing DABG calls.
#' @param thresholds the P-value thresholds to use
#' @param col the colours, 1 per threshold
#' @param grid logical: overlay a grid?
#' @param plot logical: create the plot?
#' @param legend logical: add a legend to the plot, regarding the colour
#'   thresholds?
#' @return none.
#' @author Mark Cowley
#' @export
#' @importFrom mjcgraphics colour.step axis.percentiles hgrid
hist_DABG <- function(dabg, 
					  thresholds=c(1e-05, 1e-04, 1e-03, 0.01, 0.1, 1),
					  col=colour.step("green", "red", steps=length(thresholds)), 
					  grid=TRUE, plot=TRUE, legend=TRUE) {
	counts <- t(as.data.frame(lapply(thresholds, function(thresh) apply(dabg<thresh, 2, sum)), stringsAsFactors=FALSE))
	rownames(counts) <- thresholds
	colnames(counts) <- colnames(dabg)
	# counts

	dCounts <- counts
	for(i in length(thresholds):2) {
		dCounts[i,] <- dCounts[i,] - dCounts[(i-1),]
	}
	dCounts
	
	if( plot ) {
		opar <- par(no.readonly=TRUE)
		
		# prop <- dCounts/nrow(dabg) * 100
		par(las=2, mgp=c(4,1,0))
		barplot(dCounts, col=col, main="Distribution of DABG calls", ylab="Frequency")
		axis.percentiles(side=4, max=nrow(dabg))
		if( grid ) hgrid(col="lightgrey")
		if( legend ) legend("topright", paste(sep="", "<", thresholds), fill=col, bg="white", inset=0.02)
		par( opar )
	}
	
	invisible(dCounts)
}
