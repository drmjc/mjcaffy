#' histogram of array ranks
#' 
#' plot each array as a histogram, to determine if the array has an
#' unusual proportion of measurements that are consistently the highest (or
#' lowest)
#' 
#' It will plot upto 36 arrays in one panel. Any more than this would be too
#' crowded.
#' For each gene in 'x' (ie each row), determine the ranks across each array.
#' Theoretically,
#' all measurements on the array should be independent, so each array should
#' have about the
#' same number of genes with rank 1, or rank 2 or 3 etc... If an array has lots
#' of genes that
#' are rank 1, or rank \code{ncol(x)} then there is a problem with the array /
#' normalisation
#' 
#' @param x a matrix or data.frame containing the array data
#' @param cols which columns to display? If NULL, use all the columns in x.
#' If \code{x} is provided with N columns and only some are these are specified
#' by \code{cols}, then the ranks are still calculated on all the cols in \code{x}
#' @param abline plot a horizontal line with the expected number of genes per
#' rank?
#' @param main either a single header, or 1 header per column. If NULL, then
#' colnames will be used.
#' @param xlab obvious
#' @param auto.mfrow work out what par(mfrow=c(_,_)) should be used? see help
#' for auto.mfrow.R
#' @param \dots arguments passed to \code{\link[graphics]{hist}}
#' 
#' @return invisibly returns the ranks of x - same dimensions as x.
#' 
#' @author Mark Cowley, 3 June 2005
#' @export
hist_array_ranks <- function(x, cols=NULL, abline=T, main=NULL, xlab="rank (low to high)", auto.mfrow=T, ...) {

	if(is.null(cols))
		cols <- c( 1:ncol(x) )

	stopifnot(length(cols) <= 36)

	ranks <- t(apply(x,1,rank))

	if( auto.mfrow )
		auto.mfrow(length(cols), T) ## see mjc/R/auto.mfrow.R

	if( is.null(main) )
		main <- colnames(x)
	else if( length(main)==1 & length(cols) > 1 )
		main <- paste(main, "- column", cols)

	for(i in 1:length(cols)) {
		hist(ranks[, cols[i]], breaks=c(0:max(ranks[1,]))+0.5, main=main[i], xlab=xlab, ...)
		if(abline) abline(h=nrow(x)/ncol(x), lty=2, col=2)
	}

	if( auto.mfrow )
		auto.mfrow(length(cols), F)

	invisible(ranks)
}
