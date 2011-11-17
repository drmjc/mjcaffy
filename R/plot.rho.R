## Function to work out the correlation of N slides using Spearmans rho, then plot these
## correlations as a heatmap, and also the average correlation of each array,
##
## Parameters:
##     x: The expression data; or the correlation matrix (which must thus be square with a diagonal of 1's)
##     name: Used for naming the plots appropriately
##     do.mfrow: T means call par(mfrow=c(1, 2)); F means don't set the graphical layout
##     col: the colour to use. It will be recycled if necessary
##
## Value:
##     invisibly returns the average correlation of each array.
##
## Mark Cowley, 9 November 2005
##
plot.rho <- function(x, do.mfrow=T, name=paste(ncol(x), "arrays"), col=1 ) {
    ## does col need recycling?
    col <- recycle( col )

    ## is x the expression data or a correlation matrix
    if( is.cor.matrix(x) )
        rho <- x
    else
        rho <- rho(x)

    rho.av <- rowMeans(rho)

    if( do.mfrow )
        par(mfrow=c(1,2))

    image.table(rho, main=p("rho correlations of ", name), xlab="array index", ylab="array index")
    plot(rho.av, main=p("average rho correlations of ", name), ylab="rho", xlab="array index")

    invisible( rho.av )
}
