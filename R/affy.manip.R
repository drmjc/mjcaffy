get.cell.coordinates <- function(ProbeSetID) {
    stopifnot( length(ProbeSetID) == 1 )

    if( is.numeric(ProbeSetID) )
        ProbeSetID <- geneNames(celdata)[ProbeSetID]

    system(paste("~/data/affymetrix/SpikeIn/cdf2cell.idx.sh", ProbeSetID, "> /tmp/psid.idx"))
    res <- read.delim("/tmp/psid.idx", as.is=T, header=F)
    colnames(res) <- c("X", "Y")
    for( i in seq(1,nrow(res), 2) ) {
        if(res[i+1,2] < res[i, 2]) {
            ymin <- res[i+1, 2]
            res[i+1,2] <- res[i,2]
            res[i,2] <- ymin
        }
    }
    return(res)
}

get.cell.row.idx <- function(coordinates, NROW=712, NCOL=712) {
    if(is.data.frame(coordinates) & nrow(coordinates) > 1) {
        res <- rep(0, length=nrow(coordinates))
        for(i in 1:length(res))
            res[i] <- get.cell.row.idx( coordinates[i,], NROW, NCOL )
        return(res)
    }

    return( as.numeric(coordinates[2] * NROW + coordinates[1] + 1) ) ## +1 since 0,0 is on row 1
}

get.pm.idx <- function(coordinates, NROW=712, NCOL=712) {
    idx <- get.cell.row.idx(coordinates, NROW, NCOL)
    return( idx[seq(1,length(idx), 2)] )
}

get.mm.idx <- function(coordinates, NROW=712, NCOL=712) {
    idx <- get.cell.row.idx(coordinates, NROW, NCOL)
    return( idx[seq(2,length(idx), 2)] )
}

## === interleave(pm(), mm())
get.raw <- function(ProbeSetID, array.idx) {
    coords <- get.cell.coordinates(ProbeSetID)
    idx <- get.cell.row.idx(coords)
    return(affydata@exprs[idx, array.idx])
}

## === pm()
get.raw.pm <- function(ProbeSetID, array.idx) {
    raw <- get.raw(ProbeSetID, array.idx)
    return(raw[seq(1,length(raw),2)])
}

##  === mm()
get.raw.mm <- function(ProbeSetID, array.idx) {
    raw <- get.raw(ProbeSetID, array.idx)
    return(raw[seq(2,length(raw),2)])
}



## modify.CELL.data <- function(ProbeSetID, array.idx, PM=log2(1.1), MM=0, affydata) {
##     if( length(array.idx) > 1 ) { ## the recurse through the array indices
##         for(i in 1:length(array.idx)) {
##             affydata <- modify.CELL.data(ProbeSetID, array.idx[i], PM, MM, affydata)
##         }
##         return(affydata)
##     }
## #######################
##     coords <- NULL
##     for( i in 1:length(ProbeSetID) )
##         coords <- rbind( coords, get.cell.coordinates( ProbeSetID[i] ) )
##
##     pm.idx <- get.pm.idx(coords)
##     mm.idx <- get.mm.idx(coords)
##
##     for(i in 1:length(pm.idx)) {
##         if( PM != 0 )
##             affydata@exprs[pm.idx[i], array.idx] <- 2 ^ (log2(affydata@exprs[pm.idx[i], array.idx]) + PM)
##         if( MM != 0 )
##             affydata@exprs[mm.idx[i], array.idx] <- 2 ^ (log2(affydata@exprs[mm.idx[i], array.idx]) + MM)
##     }
##
##     return( affydata )
## }



modify.CELL.all <- function(file.name, out.file.name, genes=NULL, idx=NULL, PM=1.1, MM=1.0) {
    MEAN <- 3
    ## PM on even rows seq(0, x, 2), MM on odd rows seq(0, x, 2)
    INTENSITY.idx <- as.numeric(gsub(":.*", "", system("grep \"^\\[INTENSITY\" -n CEL2/Expt3_R1.CEL", intern=T)))
    MASKS.idx <- as.numeric(gsub(":.*", "", system("grep \"^\\[MASKS\" -n CEL2/Expt3_R1.CEL", intern=T)))

    x <- read.delim(file.name, skip=INTENSITY.idx+1, as.is=T, nrow=MASKS.idx-2-INTENSITY.idx-2,
                    blank.lines.skip=F)
    ## MODIFY THE MEAN DATA BY THE FACTORS IN PM AND MM (1.1 means add 10%)
    if( is.null(genes) & is.null(idx) )
        x[,3] <- x[,3] * rep(c(PM, MM), nrow(x)/2)
    else if( !is.null(idx) ) {
        x[idx, 3] <- x[idx, 3] * rep(c(PM, MM), length(idx)/2)
    }
    else if( !is.null(genes) ){
        ## get coordinates of the specified genes.
        coords <- NULL
        for(i in 1:length(genes))
            coords <- rbind( coords, get.cell.coordinates( genes[i] ) )
        stopifnot(even(length(coords)))

        ## work out which rows of the CEL file those coords map to
        idx <- get.cell.row.idx(coords)

        x[idx, 3] <- x[idx, 3] * rep(c(PM, MM), length(idx)/2)
    }

    f <- file(out.file.name, "w")

    hdr <- scan.text(file.name, nlines=INTENSITY.idx+2)
    write(hdr, f)

    write.delim(x,f, col.names=F)

    tail <- scan.text(file.name, skip=MASKS.idx - 1)
    write.delim(tail, f)
    close(f)
}
