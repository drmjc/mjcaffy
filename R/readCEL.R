#' Function to import an Affymetrix CEL file
#' 
#' The file can be unzipped, or gzipped.
#' Since we are only interested in the [INTENSITY] section, we need to find
#' where to start
#' and stop reading the file. If the type is one of the known array types,
#' please specify
#' this in the 'type' parameter. If left as NULL, then each file will be
#' searched for when
#' to start and stop reading which takes longer.
#' 
#' @param x The file name. Can be a vector of filenames.
#' @param type MG_U74Av2 | HG_Focus | RAE230A | NULL
#' @param skip The number of lines to skip (if known)
#' @param numlines The number of lines in the [INTENSITY] section to read (if
#' known)
#' 
#' @return if x is only 1 file, then a \code{data.frame} is returned with 5 columns:
#' "X", "Y", "mean", "sd", "npixels"
#' if length(x) > 1 then a \code{list} of \code{data.frames} is returned.
#' @author Mark Cowley, 28 Oct 2005
#' @export
readCEL <- function(x, type=NULL, skip=24, numlines=NULL) {

    if( length(x) > 1 ) {
        res <- list()
        for(i in 1:length(x))
            res[[i]] <- readCEL(x[i])
        return(res)
    }
    else {
        ## File can be gzipped
        if( length( grep(".gz", x, ignore.case=T) ) > 0 ) {
            tmpname <- sub(".GZ", ".gz", paste(tempdir(), "/", basename(x), sep=""))
            file.copy(x, tmpname )
            file.gunzip( tmpname )
            cel <- readCEL( sub(".gz", "", tmpname) ) ## recursion
            unlink( sub(".gz", "", tmpname) )
        }
        else {
            ## OR unzipped
            if( is.null(skip) || is.null(numlines) ) {
                if( length(grep("Focus", type[1], ignore.case=T)) > 0 ) {
                    numlines <- 200704
                    skip <- 24
                }
                else if( length(grep("U74Av2", type[1], ignore.case=T)) > 0 ) {
                    numlines <- 409600
                    skip <- 24
                }
                else if( length(grep("RAE230A", type[1], ignore.case=T)) > 0 ) {
                    numlines <- 362404
                    skip <- 24
                }
                else {
                    ## Find the skip, numlines info
                    tmp <- scan.text(x, n=50)
                    skip <- grep("CellHeader=X", tmp)
                    numlines <- as.numeric(sub("NumberCells=", "", grep("NumberCells", tmp, value=T)))
                }
            }
            cel <- read.delim( x, as.is=T, skip=skip, header=F, nrows=numlines )
            colnames(cel) <- c("X", "Y", "mean", "sd", "npixels")
        }

        return( cel )
    }
}
