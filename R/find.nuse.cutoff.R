## Function to take a affyPLM object obtained from fitting fitPLM, and perform NUSE
## plots, identifying outlier based on various measure based on medians and IQR of
## the normalised-unscaled-standard-errors that NUSE calculates
##
## Parameters:
##     x: an rmaPLM object which has at least the standard errors stored
##     file: a pdf file to plot the nuse plots to (NULL means plot to dev.current())
##     nplots: If x has many samples, it may be useful to split the it into nplots
##     outlier.col: The colour to use for the outlying arrays
##     type: "plot" means do the plots and invisibly return the outlier arrays;
##           "outliers" means don't do the plots and return the outlier arrays.
##
## Value:
##     A list where each element is a numerical vector of indices into the columns/samples of x
##     that are 'outlier arrays'. The names of each element in the list indicate the threshold
##     test that was applied.
##
## Mark Cowley, 8 Nov 2005
##
find.nuse.cutoff <- function(x, file=NULL, nplots=1, outlier.col="red", type=c("plot", "outliers")) {
    if( type[1] == "plot" && !is.null(file) )
        pdf.A4(file)

    res <- list()

    nuse <- NUSE(x, type="stat") ## [1,] is the median, and [2,] is the IQR.

    if( nplots > 1 )
        par(mfrow=c(nplots,1), mar=c(1,4,1,2))
    par(lwd=0.7, cex.axis=0.7, las=2)

    col=rep("white", ncol(nuse))
    col[which(nuse[1,] > (mean(nuse[1,]) + sd(nuse[1,])))] <- outlier.col
    res[[length(res) + 1]] <- which(nuse[1,] > (mean(nuse[1,]) + sd(nuse[1,])))

    if( type[1] == "plot" ) {
        nuse(x, col=col, main="Threshold median above average median + 1sd of median", nplots=nplots)
        abline(h=mean(nuse[1,]) + sd(nuse[1,]), lty=3, col="purple")
        legend.TR(sum(col==outlier.col), bty="n", fill=outlier.col)
    }

    col=rep("white", ncol(nuse))
    col[which(nuse[1,] > (mean(nuse[1,])))] <- outlier.col
    res[[length(res) + 1]] <- which(nuse[1,] > (mean(nuse[1,])))

    if( type[1] == "plot" ) {
        nuse(x, col=col, main="Threshold median above average median", nplots=nplots)
        abline(h=mean(nuse[1,]), lty=3, col="purple")
        legend.TR(sum(col==outlier.col), bty="n", fill=outlier.col)
    }

    col=rep("white", ncol(nuse))
    col[which(nuse[1,] > (median(nuse[1,])))] <- outlier.col
    res[[length(res) + 1]] <- which(nuse[1,] > (median(nuse[1,])))

    if( type[1] == "plot" ) {
        nuse(x, col=col, main="Threshold median above median of the medians", nplots=nplots)
        abline(h=median(nuse[1,]), lty=3, col="purple")
        legend.TR(sum(col==outlier.col), bty="n", fill=outlier.col)
    }

    ## use IQR to determine which are the bad arrays
    col=rep("white", ncol(nuse))
    col[which(nuse[2,] > (mean(nuse[2,])))] <- outlier.col
    res[[length(res) + 1]] <- which(nuse[2,] > (mean(nuse[2,])))

    if( type[1] == "plot" ) {
        nuse(x, col=col, main="Threshold IQR above average IQR", nplots=nplots)
        legend.TR(sum(col==outlier.col), bty="n", fill=outlier.col)
    }

    col=rep("white", ncol(nuse))
    col[which(nuse[2,] > (median(nuse[2,])))] <- outlier.col
    res[[length(res) + 1]] <- which(nuse[2,] > (median(nuse[2,])))

    if( type[1] == "plot" ) {
        nuse(x, col=col, main="Threshold IQR above median IQR", nplots=nplots)
        legend.TR(sum(col==outlier.col), bty="n", fill=outlier.col)
    }

    if( type[1] == "plot" && !is.null(file) )
        dev.off()

    names(res) <- c("median > mean(medians) + sd(medians)",
                    "median > mean(medians)",
                    "median > median(medians)",
                    "IQR > mean(IQR)",
                    "IQR > median(IQR)")

    if( type[1] != "plot" )
        return(res)
    else
        invisible(res)
}
