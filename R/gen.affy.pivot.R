#' gererate Affymetrix pivot table
#' 
#' Perform RMA and either DABG or MAS5 calls depending on the Array version (3'
#' IVT or WT - whole transcript). Note, an antigenomic.bgp file MUST exist in
#' the lib.dir IF performing DABG P-values. make these files using apt-dump-pgf
#' specifying the pgf/clf and a probeset_id file. Get the antigenomic
#' probeset_id's from any of the other annotation files (the probe.tab, or the
#' annot.csv files are good places to start looking in).
#' 
#' @param cel.files a character vector of files names for each CEL file
#' @param cel.file.dir dir containing CEL files.
#' @param lib.dir (only for newer WT arrays): The location of the library files
#'   (pgf/clf/bgp). This can be either the top level library dir (eg ~/libs")
#'   or ~/libs/HuGene-1_0-st-v1). If the former, then a dir with the exact name
#'   returned by whatcdf("your.file1.CEL") must exist; i.e., a dir named
#'   HuGene-1_0-st-v1 must exist within the lib.dir.
#' @author Mark Cowley, 29/1/08
#' @export
#' @importFrom affy whatcdf justRMA read.affybatch mas5calls 
gen.affy.pivot <- function(cel.files=NULL, cel.file.dir=NULL, lib.dir="/Users/marcow/data/Microarray Libraries") {
	
	if( is.null(cel.files) ) {
		cat("Choose any CEL file from the directory that contains all CEL files.")
		file <- file.choose()
    # 		cel.file.dir <- dirname(file)
		cel.files <- dir(cel.file.dir, pattern="CEL$", full.names=TRUE)
		cat(cel.files)
	}
	else {
    # 		cel.file.dir <- dirname( cel.files[1] )
    # 		cel.files <- basename( cel.files )
	}
	
	#
	# 3' IVT or WT?
	#
	cdf <- whatcdf( cel.files[1] )
	wt.array <- grepl("gene|exon", cdf, ignore.case=TRUE)
	cat( paste("This is the:", cdf, "array, thus it is a", c("3' IVT array", "Whole Transcript (WT) Array\n")[wt.array+1]) )

	cat("Running RMA\n")
	rma <- justRMA( filenames=basename(cel.files), celfile.path=dirname(cel.files) )
	#cdf <- rma@annotation
	
	if( wt.array ) {
		# DABG summary
		cat("Running DABG\n")
		calls <- run.dabg(cel.files, lib.dir, interactive=FALSE)
	}
	else {
		cat("Running mas5calls\n")
		cel <- read.affybatch(filenames=cel.files)
		calls <- exprs(mas5calls(cel))
	}
	rma <- round(2^exprs(rma), 2)
	rows <- intersect(rownames(rma), rownames(calls))
	res <- as.data.frame(interleave.columns(rma[rows,], calls[rows,]), stringsAsFactors=FALSE)
	if( !wt.array ) {
		colclasses(res) <- c("numeric", "character")
	}
	colnames(res) <- paste(sep="", rep(sub("\\.CEL", "", colnames(rma)), each=2), c("_SIGNAL", "_DETECTION"))
	
	return( res )
}


#' run Affymetrix DABG algorithm
#' 
#' Perform DABG method to determine P-values for each ProbeSet x array
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
		# 		stop(p("lib.dir incorrectly specified\n",
		# 			   "Either specify the dir that contains the pgf/clf/bgp files,",
		# 			   "or the lib.dir that contains the ",
		# 			   shQuote(cdf), " directory"))
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

