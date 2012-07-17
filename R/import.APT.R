#' Import a file produced by APT.
#' 
#' Import a file produced by Affymetrix Power Tools (APT), (NOT debian's apt-get.org!)
#' @param file the file name
#' @param keep.first.column logical: the first column can contain (IMO) useless row ID's. 
#'    if \code{TRUE}, then include them; if \code{FALSE}, exclude this first column.
#' @param check.names logical: check the column names during the \code{\link{read.csv}},
#'    or allow them to not strictly conform to \R's naming standards (\code{FALSE})
#' @param \dots additional arguments passed to \code{\link{read.csv}}, or \code{\link{read.delim}}
#' @return a \code{data.frame} with row and col names, containing data from an APT
#'   analysis.
#' @author Mark Cowley, 16 April 2008
#' @examples
#' \dontrun{
#' rma <- import.APT("rma.summary.txt")
#' }
#' @export
import.APT <- function(file, keep.first.column=FALSE, check.names=FALSE, ...) {

	tryCatch({
		suppressWarnings({
			header <- readLines(file, 100)
			skip <- max(grep("^#", header))
		})

		if( grepl("csv$", file, ignore.case=TRUE) )
			aptfile <- read.csv(file, skip=skip, stringsAsFactors=FALSE, check.names=check.names, ...)
		else
			aptfile <- read.delim(file, skip=skip, stringsAsFactors=FALSE, check.names=check.names, ...)
	},
	error=function(e) {
		stop(simpleError("Couldn't read your file - are you sure it's a text-based Affymetrix csv file?\n"))
	})

	aptfile[aptfile=="---"] <- NA
	
	# make first column become the rownames
	rownames(aptfile) <- make.unique(as.character(aptfile[, 1]), sep=".")
	if( !keep.first.column )
		aptfile <- aptfile[, -1]
	
	return( aptfile )	
}

#' Import DABG results.
#' Import a \dQuote{dabg.summary.txt} file from running the Affymetrix
#' Detected Above Background method (DABG), skipping the header
#' @param file the path to the DABG file
#' @param \dots arguments passed to \code{\link{import.APT}}
#' @return a \code{data.frame} containing DABG P-values.
#' @seealso \code{\link{dabg2calls}} \code{\link{import.APT}}
#' @author Mark Cowley, 11/3/08
#' @export
import.dabg <- function(file="dabg.summary.txt", ...) {
	dabg <- import.APT(file, ...)
	return( dabg )	
}
