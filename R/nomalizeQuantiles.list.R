## Function to qnorm a list of matrices
## It joins the list elements into one large data.frame, then
## quantile normalizes this large data.frame.
## This then gets deconstructed back into a list with the same
## properties as the supplied list
##
## Mark Cowley, 23 Nov 2005
##
normalizeQuantiles.list <- function(x) {
    require(limma)

    tmp <- normalizeQuantiles( as.data.frame(x) )

    res <- list()
    for(i in 1:length(x)) {
        res[[i]] <- tmp[, 1:ncol(x[[i]])]
        tmp <- tmp[, -c(1:ncol(x[[i]]))]
        colnames(res[[i]]) <- colnames(x[[i]])
    }
    names(res) <- names(x)

    return(res)
}
