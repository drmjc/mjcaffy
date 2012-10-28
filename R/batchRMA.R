#' RMA normalize in batches
#' 
#' Normalize a set of arrays, using RMA in batches according to the date that
#' they were processed.
#' 
#' Each array's batch can be specified via anything that uniquely identifies 
#' arrays from the same batch, eg numerical ID, a character name, or a character date
#' 
#' @param files a vector of CEL filenames
#' @param batch a vector that identifies the batches. see Details.
#' @return a \code{data.frame} of RMA normalised data;
#' @note further normalisation may be (usually is) required between batches, eg 
#' quantile norm.
#' 
#' @examples
#' \dontrun{
#' rma.raw <- batchRMA(files, c(1,1,1,2,2,2,3,3,3))
#' rma <- normalizeQuantiles(rma.raw)
#' }
#' @author Mark Cowley, 10 April 08
#' @export
#' @importFrom affy just.rma
#' @importFrom Biobase exprs
batchRMA <- function(files, batch) {
	
	batches <- unique(batch)
	res <- list()
	for(b in batches) {
		cat("RMA-ing batch ", dQuote(b), ", files: ", dQuote(tocsv(files[batch==b])), "\n")
		res[[as.character(b)]] <- just.rma(filenames=files[batch==b])
	}
	res <- lapply(res, exprs)
	res <- cbind.list(res)
	colnames(res) <- basename(files)
	res
}
