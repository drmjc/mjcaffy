## Function to import a CEL file into an R data.frame.
##
## The file can be unzipped, or gzipped.
## Since we are only interested in the [INTENSITY] section, we need to find where to start
## and stop reading the file. If the type is one of the known array types, please specify
## this in the 'type' parameter. If left as NULL, then each file will be searched for when
## to start and stop reading which takes longer.
##
## Parameters:
##     x: The file name. Can be a vector of filenames.
##     type: MG_U74Av2 | HG_Focus | RAE230A | NULL
##     skip: The number of lines to skip (if known)
##     numlines: The number of lines in the [INTENSITY] section to read (if known)
##
##   Leave type, skip and numlines as NULL, 24, NULL for their values to be auto. determined.
##
## Value:
##     if x is only 1 file, then a data.frame is returned with 5 columns:
##        c("X", "Y", "mean", "sd", "npixels")
##     if length(x) > 1 then a list of data.frames is returned.
##
## Mark Cowley, 28 Oct 2005
##
##
## A typical CEL file looks like this....
##     ...
##     [INTENSITY]
##     NumberCells=200704
##     CellHeader=X	Y	MEAN	STDV	NPIXELS
##     0	  0	964.8	233.4	 16
##     1	  0	4954.8	480.9	 16
##     ...
##     446	447	4883.3	696.9	 16
##     447	447	63.8	12.3	 16
##
##     [MASKS]
##     NumberCells=0
##     CellHeader=X	Y
##     ...
##
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
            tmpname <- sub(".GZ", ".gz", p(tempdir(), "/", basename(x)))
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
