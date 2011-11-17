cgen2affy <- function(ids) {
    map <- read.delim("~mark/data/resourcerer/mouse/MGU74Av2_vs_Compugen22k.mastermap.txt", as.is=T)
    if( length(grep("t$",ids)) > 0 ) { ## Then ids are affy ids's -- get the genbank ids's
        res <- mgrep(ids, map[,2], nomatch=NA)
        for(i in 1:length(res)) {
            res[[i]] <- unique( map[res[[i]], 1] )
        }
    }
    else {
        res <- mgrep(ids, map[,1], nomatch=NA)
        for(i in 1:length(res)) {
            if( !is.na(res[[i]][1]) )
                res[[i]] <- unique( map[res[[i]], 2] )
        }
    }
    return( res )
}
