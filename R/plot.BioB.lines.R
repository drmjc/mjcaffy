## Plot a vertical line where the average BioB threshold is.
## The word "BioB" can be added if text.label is TRUE, and it's position relative to the line is specified by
##  the first value in label.pos.
##
## Parameters:
##    rma: use the RMA BioB threshold, or ir FALSE, use the mas5 threshold
##    text.label: add the word "BioB"??
##    label.pos: where to print the word "BioB"? The first value is chosen
##               T implies at the top of the plot (scale % from the edge)
##               B implies at the bottom of the plot (scale % from the edge)
##               L implies to left of the line, R is to the right of the line;
##    scale: used to define where the word is eg 0.05 means 5% from top (or bottom)
##    col: the colour of the line and word
##
## Value:
##    nill
##
## Mark Cowley, 3 June 2005
##
BioB.vline <- function(rma=T, text.label=T, label.pos=c("TL", "TR", "BL", "BR"), scale=0.05, col="purple", ...) {
	x <- 1
	if(rma) x <- mean(apply(expression$rma[grep("biob", rownames(expression$rma), T), -c(1:3)], 1, mean))
	else    x <- mean(apply(expression$mas5[grep("biob", rownames(expression$rma), T), -c(1:3)], 1, mean))

	abline(v=x, lty=2, col=col, ...)

	if( text.label ) {
		y <- 1
		if( contains("T", label.pos[1]) )
			y <- par("usr")[4] - diff(par("usr")[3:4])*scale
		else
			y <- par("usr")[3] + diff(par("usr")[3:4])*scale

		pos <- 2
		if( contains("L", label.pos[1]) )
			pos <- 2
		else
			pos <- 4

	## print(paste(x,y,pos))

		text(x, y, labels="BioB", pos=pos, cex=0.9, col=col, ...)
	}
}


## Plot a horizontal line where the average BioB threshold is.
## The word "BioB" can be added if text.label is TRUE, and it's position relative to the line is specified by
##  the first value in label.pos.
##
## Parameters:
##    rma: use the RMA BioB threshold, or ir FALSE, use the mas5 threshold
##    text.label: add the word "BioB"??
##    label.pos: where to print the word "BioB"? The first value is chosen
##               T implies above the line, B below the line;
##               L implies on left side of plot (scale % from the edge)
##               R implies on right side of plot (scale % from the edge)
##    scale: used to define where the word is eg 0.05 means 5% from left (or right)
##    col: the colour of the line and word
##
## Value:
##    nill
##
## Mark Cowley, 3 June 2005
##

BioB.hline <- function(rma=T,text.label=T, label.pos=c("TL", "BL", "TR", "BR"), scale=0.02, col="purple", ...) {
	y <- 1
	if(rma) y <- mean(apply(expression$rma[grep("biob", rownames(expression$rma), T), -c(1:3)], 1, mean))
	else    y <- mean(apply(expression$mas5[grep("biob", rownames(expression$rma), T), -c(1:3)], 1, mean))

	abline(h=y, lty=2, col=col, ...)

	if( text.label ) {
		x <- 1
		if( contains("L", label.pos[1]) )
			x <- par("usr")[1] + diff(par("usr")[1:2])*scale
		else
			x <- par("usr")[2] - diff(par("usr")[1:2])*scale

		pos <- 2
		if( contains("T", label.pos[1]) )
			pos <- 3
		else
			pos <- 1

	## print(paste(x,y,pos))

		text(x, y, labels="BioB", pos=pos, cex=0.9, col=col, ...)
	}
}

