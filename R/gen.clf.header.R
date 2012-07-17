#%chip_type=HuEx-1_0-st-v2
#%chip_type=HuEx-1_0-st-v1
#%chip_type=HuEx-1_0-st-ta1
#%lib_set_name=HuEx-1_0-st
#%lib_set_version=r2
#%create_date=Tue Sep 19 15:18:05 PDT 2006
#%guid=0000008635-1158704285-0732263232-1857033251-0689718480
#%clf_format_version=1.0
#%rows=2560
#%cols=2560
#%sequential=1
#%order=col_major
#%header0=probe_id	x	y

# Change the lib_set info for a CLF file.
#  NB, this differs in function to most of the other gen.xxx.header style methods
# because we keep the same CLF that belongs to the chip type of the CEL files,
# whereas in the pgf/mps files, we are changing the header of the EXON files to
# point towards data from a GENE CEL file.
#
# Mark Cowley, 7/1/08
# #
# change.clf <- function(in, out, lib_set_name="HuExonLite-1_0-st-v1", lib_set_version="v1")
# 	rename.output <- FALSE
# 	if( in == out ) {
# 		warning("This will destroy the original CLF file\n")
# 		out <- tempfile()
# 		rename.output <- TRUE
# 	}
# 	file.copy(in, out)
# 	
# 	if( rename.output )
# 		file.move(out, in)
# }

#' gen.clf.header
#' 
#' Generate a CLF file header.\cr
#' CLF is the chip layout file from Affymetrix
#'
#' @param chip_type A vector supported Affymetrix array type(s)
#' @param lib_set_name The custom lib_set_name within the header
#' @param version The custom lib_set_version
#' @param clf.header.template A named vector from reading in a header.
#' @return a vector corresponding to a CLF file header
#' @author Mark Cowley, 2011-10-31
#' @export
gen.clf.header <- function(chip_type=c("HuGene-1_0-st-v1"), lib_set_name="HuExonLite-1_0-st-v1", version="v1", clf.header.template=NULL) {
	if( is.null(clf.header.template) )
		res <- NULL
	for(ct in chip_type)
		res <- c(res, paste(sep="", '#%chip_type=', ct))

	res <- c(	res,
				p('#%lib_set_name=', lib_set_name),
				p('#%lib_set_version=', version),
				p('#%create_date=', date()),
				p('#%guid=', gen.affy.guid()),
				p('#%clf_format_version=', clf.header.template$clf_format_version),
				p('#%rows=',clf.header.template$rows),
				p('#%cols=',clf.header.template$cols),
				p('#%sequential=',clf.header.template$sequential),
				p('#%order=', clf.header.template$order),
				p('#%header0=', clf.header.template$header0)
			)
	return( res )
}


#' parse out the info in the header of a CLF file.
#' 
#' @param clf.file the path to a clf file
#' @return a named list of values within a CLF file header
#' @author Mark Cowley, 7/1/08
#' @export
#' @examples
#' \dontrun{
#' parse.clf.header("~/data/Microarray Libraries/HuGene-1_0-st-v1/HuGene-1_0-st-v1.r3.clf")
#' parse.clf.header("~/data/Microarray Libraries/HuEx-1_0-st-v2/HuEx-1_0-st-v2.r2.clf")
#' }
parse.clf.header <- function(clf.file) {
	header <- readLines(clf.file, 30)
	header <- header[grep("^#%", header)]
	header <- sub("#%", "", header)
	keys <- sub("=.*", "", header)
	vals <- sub("^.*=", "", header)
	map <- cbind(keys=keys, vals=vals)
	res <- map2list(map)
    # 	res <- as.list(vals)
    # 	names(res) <- keys
	
	return( res )
}
