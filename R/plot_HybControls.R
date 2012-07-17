
plot.HybControls <- function(x, col=NULL, lty=NULL, legend=TRUE, ...) {
    names <- c("BioB", "BioC", "BioD", "cre")
    Econc <- c(1.5,5,25,100)

	if( is.null(col) )
		col <- 1:8
	col <- recycle(col, ncol(x))
	if( is.null(lty) )
		lty <- 1
	lty <- recycle(lty, ncol(x))

    tmp <- matrix(NA, ncol(x), length(names))
    colnames(tmp) <- names
    rownames(tmp) <- colnames(x)

    for(i in 1:length(names))
        tmp[,i] <- apply(x[grep(names[i], rownames(x), TRUE),], 2, mean)

    plot(Econc, tmp[1,], type="l", ylim=c(0, max(x, na.rm=TRUE)),
         xlab="Expected concentration (pM)",
         ylab="Expression level", main="Hybridisation Controls plot", xaxt="n",
		col=col[1], lty=lty[1], ...)
    axis(side=1, at=Econc, labels=Econc)
    text(Econc, c(0,0.5,0,0), labels=names, cex=0.8)
    for(i in 2:nrow(tmp))
        lines(Econc, tmp[i,], col=col[i], lty=lty[i])

	if( legend ) {
		legend("bottomright", colnames(x), 
			col=col, lty=lty, inset=0.05, ncol=floor(ncol(x)/25)+1, cex=0.5)
	}

    invisible(tmp)
}
## tmp2 <- plot.HybControls(expression.raw$rma)

