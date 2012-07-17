#' run Affymetrix DABG algorithm
#' 
#' Perform DABG method to determine P-values for each ProbeSet x array
#' 
#' @section TODO:
#' use \code{\link{Sys.which}}
#' 
#' @param cel.files a character vector of cel file names
#' @param lib.dir the location of the PGF/CLF files; or set \code{NULL} & specify \code{interactive=TRUE}
#' @param interactive logical: if \code{TRUE}, the ask user for the CEL file dir, PGF, CLF, BGP files
#' 
#' @return the DABG calls as a \code{data.frame}
#' 
#' @author Mark Cowley, 29/1/08
#' @export
run.dabg <- function(cel.files=NULL, lib.dir=NULL, interactive=FALSE) {
	
	apt <- system("which apt-probeset-summarize", intern=TRUE)
	if( length(apt) == 0 ) {
		msg <- c("Couldn't find your installation of apt-probeset-summarize",
				 "try: system(\"echo $PATH\")",
				 "and make sure that the apt-probeset-summarise is in one of those locations")
		stop( msg )
	}
	
	if( interactive ) {
		cat("Choose any CEL file from the directory that contains all CEL files\n")
		cel.file.dir <- dirname(file.choose())
		cel.files <- dir(cel.file.dir, pattern="CEL*", full.names=TRUE)
		cat("Choose the PGF file\n")
		pgf <- file.choose()
		cat("Choose the CLF file\n")
		clf <- file.choose()
		cat("Choose the antigenomic BGP file\n")
		bgp <- file.choose()
		cat("Choose the probeset kill list file\n")
		kill.list <- file.choose()
	}
	else {
		# # is this the top level lib.dir or the level where the PGF/CLF files are?
		# if( length(dir(lib.dir, pattern="pgf$")) == 0 ) {
		# 	cdf <- cdfName( cel.files[1] )
		# 	lib.dir <- file.path(lib.dir, cdf)
		# 	if( length(dir(lib.dir, pattern="pgf$")) == 0 ) {
		# 		stop(paste("lib.dir incorrectly specified\n",
		# 			   "Either specify the dir that contains the pgf/clf/bgp files,",
		# 			   "or the lib.dir that contains the ",
		# 			   shQuote(cdf), " directory", sep=""))
		# 	}
		# }
		
		cdf <- cdfName( cel.files[1] )
		pgf <- dir(lib.dir, pattern=".*\\.pgf$", full.names=TRUE)[1]
		clf <- dir(lib.dir, pattern=".*\\.clf$", full.names=TRUE)[1]
		bgp <- dir(lib.dir, pattern=".*antigenomic.*\\.bgp$", full.names=TRUE)[1]
		kill.list <- dir(lib.dir, pattern=".*.kill.list.*", full.names=TRUE)[1]
	}

	outdir <- tempdir()
	cels <- paste(shQuote(cel.files), collapse=" ")
	if( file.exists(kill.list) )
		kill.list <- paste("--kill-list", shQuote(kill.list))
	else
		kill.list <- ""
	cmd <- sprintf("%s -a dabg -p %s -c %s -b %s %s -o %s %s",
		shQuote(apt), shQuote(pgf), shQuote(clf), shQuote(bgp), kill.list, shQuote(outdir), shQuote(cels))

	# cmd <- paste(apt, "-a dabg -p",shQuote(pgf),"-c",shQuote(clf),"-b",shQuote(bgp),"-o",shQuote(outdir), cels)
	system(cmd, wait=TRUE)
	
	if( !file.exists(file.path(outdir, "dabg.summary.txt")) )
		stop("DABG file: dabg.summary.txt does not appear to exist")
	calls <- import.APT(file.path(outdir, "dabg.summary.txt"))
	# calls <- column2rownames(read.delim(file.path(outdir, "dabg.summary.txt"), as.is=TRUE, comment.char="#"), 1)
	unlink(outdir)
	
	return( calls )
}
