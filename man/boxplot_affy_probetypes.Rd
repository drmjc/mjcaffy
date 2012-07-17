\name{boxplot_affy_probetypes}
\alias{boxplot_affy_probetypes}
\title{Affymetrix boxplot of probe types}
\usage{
  boxplot_affy_probetypes(data, annot, merge = TRUE,
    main = "Expression levels for different types of probes - RMA",
    ...)
}
\arguments{
  \item{data}{a data.frame of expression values, with
  rownames.}

  \item{annot}{a data.frame of annotations, containing the
  column "category", and with rownames.}

  \item{merge}{if TRUE, then one plot is made, containing
  the rowMeans of genes in data. If false, then n plots are
  made, for each of the n columns in data.}

  \item{main}{plot title}

  \item{\dots}{other arguments passed to plot}
}
\description{
  Perform a boxplot of expression levels for all genes on
  the array, with one box per probe type.
}
\note{
  For instance, on the new HuGene_1.0st chip there are:
  28869 main probes 2904 "neg_contr" (normgene->intron)
  1195 "pos_contr" (normgene->exon) 57 control->affx 45
  control->bgp->antigenomic 227 rescue->FLmRNA->unmapped
}
\author{
  Mark Cowley, 2011-07-14
}

