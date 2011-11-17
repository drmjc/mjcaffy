## require(methods)
## require(graphics)
## require(affyPLM)
## ## plot NUSE by splitting up the number of arrays into nplots.
## setMethod("boxplot",signature(x="PLMset"),
##           function(x,type=c("NUSE","weights","residuals", "RLE"),range=0,
##                    nplots=1,col="white",...) {
##
##
##     compute.nuse <- function(which) {
##         nuse <- apply(x@weights[which,],2,sum)
##         1/sqrt(nuse)
##     }
##
##
##     type <- match.arg(type)
##     model <- x@model.description$modelsettings$model
##     if (type == "NUSE") {
##         if (x@model.description$R.model$which.parameter.types[3] == 1 & x@model.description$R.model$which.parameter.types[1] == 0 ) {
##             grp.rma.se1.median <- apply(se(x), 1,median,na.rm=TRUE)
##             grp.rma.rel.se1.mtx <- sweep(se(x),1,grp.rma.se1.median,FUN='/')
##             if(nplots==1)
##                 boxplot(data.frame(grp.rma.rel.se1.mtx),range=range,col=col,...)
##             else {
##                 ids <- split.into.batches(1:ncol(grp.rma.rel.se1.mtx), N=nplots)
##                 if(length(col) != ncol(grp.rma.rel.se1.mtx))
##                     col <- rep(col, ncol(grp.rma.rel.se1.mtx))[1:ncol(grp.rma.rel.se1.mtx)]
##                 for(i in 1:length(ids)) {
##                     boxplot(data.frame(grp.rma.rel.se1.mtx[, ids[[i]] ]), range=range, col=col[ ids[[i]] ], ...)
##                     abline(h=1)
##                 }
##             }
##         }
##         else {
##             # not the default model try constructing them using weights.
##             which <-indexProbesProcessed(x)
##             ses <- matrix(0,length(which) ,4)
##
##             for (i in 1:length(which))
##                 ses[i,] <- compute.nuse(which[[i]])
##
##
##             grp.rma.se1.median <- apply(ses, 1,median)
##             grp.rma.rel.se1.mtx <- sweep(ses,1,grp.rma.se1.median,FUN='/')
##             boxplot(data.frame(grp.rma.rel.se1.mtx),range=range,...)
##         }
##     }
##     else if (type == "weights") {
##         ow <- options("warn")
##         options(warn=-1)
##         if (x@model.description$R.model$response.variable == -1) {
##             boxplot(data.frame(x@weights[[2]]),...)
##         }
##         else if (x@model.description$R.model$response.variable == 1) {
##             boxplot(data.frame(x@weights[[1]]),...)
##         }
##         else {
##             boxplot(data.frame(rbind(x@weights[[1]],x@weights[[2]])),...)
##         }
##         options(ow)
##     }
##     else if (type == "residuals") {
##         ow <- options("warn")
##         options(warn=-1)
##         if (x@model.description$R.model$response.variable == -1) {
##             boxplot(data.frame(x@residuals[[2]]),...)
##         }
##         else if (x@model.description$R.model$response.variable == 1) {
##             boxplot(data.frame(x@wresiduals[[1]]),...)
##         }
##         else {
##             boxplot(data.frame(rbind(x@residuals[[1]],x@residuals[[2]])),...)
##         }
##         options(ow)
##     }
##     else if ( type == "RLE") {
##         RLE(x, type="plot", nplots=nplots, col=col, range=range, ...)
##     }
## })
## ## rm(ow, grp.rma.se1.median, grp.rma.rel.se1.mtx, which, ses, i, ids, type, model )
## ## rm(last.warning)
##
##
## if (!isGeneric("RLE"))
##   setGeneric("RLE",function(x,...)
##              standardGeneric("RLE"))
##
##
##
## ## New version of RLE to allow splitting the phenotypes into 'nplots' plots.
## setMethod("RLE",signature(x="PLMset"),
##             function(x,type=c("plot","values","stats"),ylim=c(-0.75,0.75),add.line=TRUE,
##                      nplots=1, col="white", range=0, ...) {
##
##
##     type <- match.arg(type)
##     model <- x@model.description$modelsettings$model
##     if (type == "values" || type=="stats") {
##         if (x@model.description$R.model$which.parameter.types[3] == 1) {
##             medianchip <- apply(coefs(x), 1, median)
##             if (type == "values") {
##                 sweep(coefs(x),1,medianchip,FUN='-')
##             }
##             else {
##                 RLE <- sweep(coefs(x),1,medianchip,FUN='-')
##                 Medians <- apply(RLE,2,median)
##                 Quantiles <- apply(RLE,2,quantile,prob=c(0.25,0.75))
##                 RLE.stats <- rbind(Medians,Quantiles[2,] - Quantiles[1,])
##                 rownames(RLE.stats) <- c("median","IQR")
##                 RLE.stats
##             }
##         }
##         else {
##             stop("It doesn't appear that a model with sample effects was used.")
##         }
##     }
##     else {
##         if(nplots==1) {
##             Mbox(x,ylim=ylim,col=col,...)
##             if (add.line)
##                 abline(0,0)
##         }
##         else { ## Changed the code from here on...
##             ids <- split.into.batches(1:ncol(coefs(x)), N=nplots)
##             if( length(col) < ncol(coefs(x)) )
##                 col <- rep(col, ncol(coefs(x)))[1:ncol(coefs(x))]
##
##             rle <- RLE(x, type="values")
##             ## recapitulate what Mbox does (calc medianchip then M then boxplot M setting the range and ylim)
##             medianchip <- apply(coefs(x), 1, median)
##             M <- sweep(coefs(x),1,medianchip,FUN='-')
##             for( i in 1:length(ids) ) {
##                 boxplot(as.data.frame(M[, ids[[i]] ]), ylim=ylim, range=range, col=col[ ids[[i]] ], ...)
##                 if (add.line)
##                     abline(0,0)
##             }
##         }
##     }
## })
