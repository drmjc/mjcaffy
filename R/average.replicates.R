## Function to summarise a table of expression measures where there might be > 1 array per strain.
## Eg if there were 100 arrays, and each array belonged to a strain, and some strains had > 1 array,
## then this function will work out the unique strains, then for each of the unique strains, look in
## strain.labels to work out which columns in expr to look in. Then it will summarise the N (>=1) arrays
## for that strain and produce a matrix of summarised expression measures.
##
## Parameters:
##   expr: a matrix or exprSet
##   strains: If colnames of expr look like:
##      strainA.1, strainA.2, strainB.1, strainB.2;
##      then you can just set strains=c("strainsA", "strainsB"), and grep
##   factor: a numeric vector indicating which columns belong to the same
##           strain. It's recommended that you also set strain to the unique
##           strain names so that the final column names are sensible. This
##           is the most reliable method:
##        eg: 9 columns for 3 replicates of 3 strains:
##        average.replicates(x, strains=LETTERS[1:3], factor=rep(1:3, each=3))
##   method: mean or median (as a character string - not the function name).
##
## Value:
##  a matrix with same nrow as expr, but one column per strain. If a strain is
##  not found in expr then it is removed from the returned results with no warning.
##
## Mark Cowley, 12 July 2005
##
average.replicates <- function(expr, strains=NULL, factor=NULL, method=c("mean", "median")) {
	stopifnot( method[1] %in% c("mean", "median") )
	stopifnot( length(strains) <= ncol(expr) )
	if( !is.null(strains) && is.null(factor) ) {
        factor <- mgrep(p(strain,"$|", strain, "[.]"), colnames(expr))
        # convert the factor list into a vector, such that all of the
        # arrays from the first strain contain idx 1. If a strain is
        # not represented, then that index will not be present in the
        # factor vector.
        factor <- rep(1:length(factor), each=sapply(factor, l))
	}

    if( class(expr) == "exprSet" )
        exprs <- expr@exprs

	res <- matrix(NA, nrow(expr), length(strains))

	for(i in 1:length(strains)) {
		strain <- strains[i]
        idx <- which(factor == i)
# 		idx <- grep(p(strain,"$|", strain, "[.]"), colnames(expr))
		if( length(idx) == 0 ) {
			res[,i] <- rep(NA, nrow(expr))
		}
		else if( length(idx) == 1 ) {
			res[,i] <- expr[,idx]
		}
		else { #}if ( length(idx) > 1 ) {
			if( method[1]=="mean" )
				res[,i] <- rowMeans(expr[,idx])
			else
				res[,i] <- rowMedians(expr[,idx])
		}
	}

    if( !is.null(strains) )
	   colnames(res) <- strains
	rownames(res) <- rownames(expr)

	for(i in ncol(res):1) {
		if(all(is.na(res[,i])))
			res <- res[,-i]
	}

	return(res)
}
