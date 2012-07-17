
plot.PolyA.Controls <- function(x, col=NULL, lty=NULL, legend=TRUE, ...) {
    names <- c("lys", "phe", "thr", "dap")
    Eratio <- c(1/100000, 1/50000, 1/25000, 1/7500)
    labels <- c("1:100,000", "1:50,000", "1:25,000", "1:7,500")

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

    plot(Eratio, tmp[1,], type="l", ylim=c(0, max(x, na.rm=TRUE)),
         xlab="Expected final concentration of controls relative to total RNA population",
         ylab="Expression level", main="Poly-A Controls plot", xaxt="n", col=col[1], lty=lty[1], ...)
    axis(side=1, at=Eratio, labels=labels)
    text(Eratio, rep(0,length(names)), labels=names, cex=0.8)

    for(i in 2:nrow(tmp))
        lines(Eratio, tmp[i,], col=col[i], lty=lty[i])

	if( legend ) {
		legend("bottomright", colnames(x), 
			col=col, lty=lty, inset=0.05, ncol=floor(ncol(x)/25)+1, cex=0.5)
	}

    invisible(tmp)
}

## plot.PolyA.Controls(expression.raw$rma)
