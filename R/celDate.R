#' datestamp of Affymetrix CEL file
#' Extract the CEL file creation date stamp from within the CEL file header.
#' 
#' @param files a vector of Affymertrix CEL files.
#' @return a vector of datestamps from the CEL file headers
#' @author Mark Cowley, 2008-07-29
#' @export
celDate <- function(files) {
	stopifnot( all(file.exists(files)) )
	files <- paste(sQuote(files), collapse=" ")
	#cmd <- paste("~/bin/celDate.sh", files)
	#grep -m1 -a '^DatHeader' "$@" | egrep -o '[0-9]{2}/[0-9]{2}/[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}'
	cmd <- paste("grep -m1 -a '^DatHeader'", files, 
				"| egrep -o '[0-9]{2}/[0-9]{2}/[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}'")
    
	dates <- system(cmd, intern=TRUE)
	dates
}
