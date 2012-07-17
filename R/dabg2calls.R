#' Convert DABG P-values into P/M/A calls, by hard-thresholding at two
#' threshold, one for "Present", and "Marginal"
#' 
#' @param dabg a data.frame of dabg calls. see ?import.APT
#' @param Pthresh the p-values at which to threshold the DABG p-values
#' @param Mthresh the p-values at which to threshold the DABG p-values
#' @return a data.frame, with same dim, and dimname as dabg, but with values of
#'   "P", "M", "A"
#' @author Mark Cowley, 2008-05-26
#' @export
dabg2calls <- function(dabg, Pthresh=1e-05, Mthresh=0.001) {
	calls <- as.data.frame(matrix("A", nrow(dabg), ncol(dabg)), stringsAsFactors=FALSE)
	dimnames(calls) <- dimnames(dabg)
	calls[dabg < Pthresh] <- "P"
	calls[dabg < Mthresh & calls != "P"] <- "M"
	
	calls
}
