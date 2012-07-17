#' Convert Affymetrix MAS5 calls into colours.
#' 
#' @param x	If \code{x} is a \code{vector} of calls, then P/M/A = green/orange/red;
#'  else if \code{x} is a \code{matrix} of calls, then allP = green, noneP = red, else orange
#' 
#' @return a character vector or colours same length as \code{length(x)} or \code{nrow(x)}
#' 
#' @author Mark Cowley, 2009-05-15
#' @export
calls2colour <- function(x) {
	green <- "#7F9A48"
	orange <- "#F79646"
	red <- "#C0504D"
	grey <- "#BEBEBE"
	
	cols <- NULL
	if( is.vector(x) ) {
		cols <- rep(grey, length(x))
		
		cols[x=="P"] <- green
		cols[x=="M"] <- orange
		cols[x=="A"] <- red
	}
	else {
		cols <- rep(grey, nrow(x))

		Pcount <- rowSums(x=="P")
		cols[Pcount==ncol(x)] <- green
		cols[Pcount>0 & Pcount<ncol(x)] <- orange
		cols[Pcount==0] <- red
	}
	cols
}
