#' Affymetrix boxplot of probe types
#' 
#' Perform a boxplot of expression levels for all genes on the array, with one
#' box per probe type.
#' 
#' @param data a data.frame of expression values, with rownames.
#' @param annot a data.frame of annotations, containing the column "category",
#'   and with rownames.
#' @param merge if TRUE, then one plot is made, containing the rowMeans of
#'   genes in data. If false, then n plots are made, for each of the n columns
#'   in data.
#' @param main plot title
#' @param \dots other arguments passed to plot
#' @note For instance, on the new HuGene_1.0st chip there are: 28869 main
#'   probes 2904 "neg_contr" (normgene->intron) 1195 "pos_contr"
#'   (normgene->exon) 57 control->affx 45 control->bgp->antigenomic 227
#'   rescue->FLmRNA->unmapped
#' @author Mark Cowley, 2011-07-14
#' @export
boxplot_affy_probetypes <- function(data, annot, merge=TRUE, main="Expression levels for different types of probes - RMA", ...) {

    tmp <- intersect(rownames(data), rownames(annot))
    data <- data[tmp,]
    annot <- annot[tmp,]

    if( merge && ncol(data) > 1 ) {
        data <- rowMeans(data)
    }

    if (merge) {
        boxplot(data ~ as.factor(annot$category), width=sqrt(ucounts(annot$category)), ...)
        title(main=main)
    }
    else {
        auto.mfrow(ncol(data))
        for(i in 1:ncol(data)) {
            boxplot_affy_probetypes(data, annot, merge=TRUE, main=paste(main, "\n", colnames(data)[i], sep=""), ...)
        }
    }
}


