#' Extract the CDF file name from a CEL file header.
#' 
#' affy provides \code{cdfName,AffyBatch}, this provides \code{cdfName,character}, where
#' a character vector of file paths can be used to obtain the CDF name embedded
#' in the CEL file header.
#' 
#' @param object either an \code{\link[affy]{AffyBatch}}, or character vector of filenames
#' 
#' @return a vector of cdf names for each CEL file.
#' 
#' @author Mark Cowley, 2008-05-26
#' 
#' @export
#' @importClassesFrom affy AffyBatch
#' @importMethodsFrom affy cdfName
#' @importFrom affyio read.celfile.header
#' 
#' @rdname cdfName-methods
#' @aliases cdfName,character-method
#' 
#' @examples
#' \dontrun{
#' files <- dir("/Volumes/Volumes/PWBC/private/projects/ShaneGrey/CEL file repsitory/CEL/", full.names=TRUE, pattern="AW.*CEL")[1:3]
#' cdfName(files)
#' }
setMethod("cdfName",
	signature=signature("character"),
	function(object) {
		filenames <- object

		res <- rep(NA, length(object))
		for( i in 1:length(object) ) {
			filename <- object[i]
			# # This snippet came from looking into affy::read.affybatch
			# try(res[i] <- .Call("ReadHeader", as.character(filename), PACKAGE = "affyio")[[1]]) # [[2]] would be the dimensions.
			res[i] <- read.celfile.header(as.character(filename))$cdfName
		}
	
		return( res )
		
	}
)
# CHANGELOG
# 2012-07-06: changed to S4 method. dropped the .Call

# cdfName <- function(object) {
# 	if( class(object) == "AffyBatch" ) {
# 		affy::cdfName(object)
# 	}
# 	else if( class(object) == "character" ) {
# 		filenames <- object
# 
# 		res <- rep(NA, length(filenames))
# 		for( i in 1:length(filenames) ) {
# 			filename <- filenames[i]
# 			# This snippet came from looking into affy::read.affybatch
# 			try(res[i] <- .Call("ReadHeader", as.character(filename), PACKAGE = "affyio")[[1]]) # [[2]] would be the dimensions.
# 		}
# 	
# 		return( res )
# 	}
# 	else {
# 		stop("Unsupported object type.\n")
# 	}
# }
