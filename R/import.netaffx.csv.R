#' import Affymetrix NetAffx csv file
#' 
#' import a file from NetAffx, designed to work with annotation
#' files that have comment rows, and are csv's.
#' Warning, this can take a long time! 41s to read the HuEx na32 transcript.csv 
#' from an SSD on a 2010 MBP, or 192s to read the HuEx na32 probeset.csv file
#' 
#' @param file the path to a NetAffx csv file
#' @param \dots other arguments passed to \code{\link{read.csv}}
#' @return a \code{data.frame}
#' @author Mark Cowley, 2008-05-27
#' @export
import.netaffx.csv <- function(file, ...) {
	header <- readLines(file, 100)
	commentlines <- grep("^#", header)
	skip <- ifelse(length(commentlines) == 0, 0, max(commentlines))

	netaffxfile <- read.csv(file, skip=skip, ...)

	# make first column become the rownames
	rownames(netaffxfile) <- make.unique(as.character(netaffxfile[, 1]), sep=".")
    # netaffxfile <- netaffxfile[, -1]
	
	return( netaffxfile )	
}
