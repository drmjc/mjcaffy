#' Import an Affymetrix MPS file
#' 
#' @param file the file name
#' @param \dots further arguments passed to read.delim (hint nlines=100)
#' @return a data.frame of 4 columns: probeset_id, transcript_cluster_id,
#'   probeset_list, probe_count
#' @author Mark Cowley, 2/1/2008
#' @export
#' @examples
#' \dontrun{
#' ex.mps <- import.mps("~/data/Microarray Libraries/HuEx-1_0-st-v2/HuEx-1_0-st-v2.r2.dt1.hg18.core.mps", nrows=100)
#' head(ex.mps)
#' }
import.mps <- function(file, ...) {
    header <- readLines(file, 30)
    skip <- max(grep("^#",header))
    mps <- read.delim(file, skip=skip, ...)
    mps$probeset_list <- trim(mps$probeset_list)
    rownames(mps) <- mps$probeset_id

    return(mps)
}

#' Convert an an Affymetrix MPS table into a named list, mapping affymtrix
#' meta-probeset_id's (the list element names) to the probeset_id's that make
#' up the meta-probeset (or the transcript cluster)
#' 
#' @param x an mps object. see import.mps
#' @return a list of N entries, corresponding to N metaprobesets, where each
#'   element is a vector of probeset_id's.
#' @author Mark Cowley, 2/1/2008
#' @export
#' @examples
#' \dontrun{
#' ex.mps2 <- mps2list(ex.mps)
#' }
mps2list <- function(x) {
    res <- strsplit(x$probeset_list, " ")
    names(res) <- x$probeset_id
    return( res )
}

#' Convert an Affymetrix MPS object into a 2column map, where each row contains
#' a key and a value; keys being the meta-probeset_id, and value being the
#' probeset_id's
#' 
#' @param x an mps object. see import.mps
#' @return a map object: data.frame with 2 columns (metaprobeset_id and
#'   probeset_id)
#' @author Mark Cowley, 2/1/2008
#' @export
#' @examples
#' \dontrun{
#' ex.mps3 <- mps2map(ex.mps)
#' head(ex.mps3)
#' }
mps2map <- function(x) {
    l <- mps2list(x)
    map <- data.frame(rep(names(l), times=sapply(l,length)), unlist(l), stringsAsFactors=FALSE)
    colnames(map) <- c("metaprobeset_id", "probeset_id")
    rownames(map) <- NULL
#    colclasses(map) <- "integer"
    return(map)
}

#' Convert an Affymetrix MPS object into a hash table (an R environment); The
#' keys being the meta-probeset_id, and value being the probeset_id's
#' 
#' @param x an mps object. see import.mps
#' @return a hash table/R environment
#' @author Mark Cowley, 2/1/2008
#' @export
mps2env <- function(x) {
    map <- mps2map(x)
    map[,1] <- as.character(map[,1])
    env <- map2env(map, 1, 2)
    return( env )
}
