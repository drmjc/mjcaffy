#' Generate a version 1.0 PGF header
#' 
#' @author Mark Cowley
#' @export
#' @examples
#' \dontrun{
#' write(gen.pgf.header(chip_type="HuGene-1_0-st-v1",
#' lib_set_name="HuExonLite-1_0-st-v1", version="v1"), f)
#' write(gen.pgf.header(chip_type="MoGene-1_0-st-v1",
#' lib_set_name="MoExonLite-1_0-st-v1", version="v1"), f)
#' write(gen.pgf.header(chip_type="RaGene-1_0-st-v1",
#' lib_set_name="RaExonLite-1_0-st-v1", version="v1"), f)
#' }
gen.pgf.header <- function(chip_type="HuGene-1_0-st-v1", lib_set_name="HuExonLite-1_0-st", version="v1") {
    res <- c( p('#%chip_type=', chip_type), 
              p('#%lib_set_name=', lib_set_name), 
              p('#%lib_set_version=', version), 
              '#%pgf_format_version=1.0', 
              paste('#%guid=', gen.affy.guid(), sep=""), 
              p('#%create_date=',system("date", intern=TRUE)), 
              '#%header0=probeset_id\ttype', 
              '#%header1=\tatom_id', 
              '#%header2=\t\tprobe_id\ttype\tgc_count\tprobe_length\tinterrogation_position\tprobe_sequence')
    res
}

#' Generate a header for an MPS file
#' 
#' todo: add mm8, mm9, rn? support.
#' 
#' @author Mark Cowley, 7/1/07
#' @export
#' @examples
#' \dontrun{
#' gen.mps.header(chip_type="HuGene-1_0-st-v1",
#' lib_set_name="HuExonLite-1_0-st-v1", version="v1", genome="hg18")
#' gen.mps.header(chip_type="MoGene-1_0-st-v1",
#' lib_set_name="MoExonLite-1_0-st-v1", version="v1", genome="mm9")
#' gen.mps.header(chip_type="RaGene-1_0-st-v1",
#' lib_set_name="RaExonLite-1_0-st-v1", version="v1", genome="rn4")
#' }
gen.mps.header <- function(chip_type="HuGene-1_0-st-v1", lib_set_name="HuExonLite-1_0-st-v1", version="v1", genome=c("hg18")) {
    res <- c( p('#%chip_type=', chip_type), 
              p('#%lib_set_name=', lib_set_name), 
              p('#%lib_set_version=', version), 
              p('#%create_date=',date()) )
    if( genome == "hg18" )
        res <- c(res, 
                 c("#%genome-species=Homo sapiens",
                    "#%genome-version=hg18",
                    "#%genome-version-ucsc=hg18",
                    "#%genome-version-ncbi=36",
                    "#%genome-version-create_date=2006 March"))
    else {
        stop("unsupported genome. please edit gen.mps.header function, inside of gen.pgf.header.R")
    }
    res <- c(res, c(paste('#%guid=', gen.affy.guid(), sep=""), 
                    "probeset_id\ttranscript_cluster_id\tprobeset_list\tprobe_count"))

    return( res )
}
