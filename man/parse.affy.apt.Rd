\name{parse.affy.apt.mrna}
\alias{parse.affy.apt.gene_assignment}
\alias{parse.affy.apt.mrna}
\title{parse Affymetrix gene info tables for ST arrays}
\usage{
  parse.affy.apt.mrna(f = NULL, annot = NULL, n = 100)

  parse.affy.apt.gene_assignment(f = NULL, annot = NULL,
    n = 100)
}
\arguments{
  \item{f}{an optional file name to the csv annotation
  file}

  \item{annot}{a pre-imported annotation data.frame. Must
  supply one of f or annot must contain at least a
  \code{probeset_id} and \code{mrna_assignment} or
  \code{gene_assignment} columns}

  \item{n}{if \code{NULL}, then process all rows of annot,
  otherwise, process first n rows.}
}
\value{
  a \code{data.frame} with 'n' rows for each row of annot,
  and many columns (see description).
}
\description{
  parse the mrna_assignment or gene_assignment column out
  of an ST transcript.csv file (an presumably a
  probeset.csv file, though this hasn't been checked).
}
\details{
  These files contain a number of different sources of
  annotation, including: RefSeq, GenBank, ENSEMBL (of
  various grades) etc... Each annotation source found in
  the annot file will be merged into 2 columns in the
  output, called eg: "ENSEMBL Transcript ID" and "ENSEMBL
  Transcript DESC", or "RefSeq ID" and "RefSeq Desc". Thus
  there will be 1+2n columns in output for the probeset_id,
  then each of the n annotation sources. When multiple
  ID/Descriptions from the same annotation source are
  identified for each gene, then they will be seperated by
  " // ".

  The nitty gritty. This needs parsing up into bits:
  annot$mrna_assignment[1] [1] "GENSCAN00000069186 //
  ENSEMBL Prediction // cdna:Genscan
  chromosome:NCBIM37:1:3018707:3044814:1 // chr1 // 100 //
  100 // 33 // 33 // 0 /// GENSCAN00000048429 // ENSEMBL
  Prediction // cdna:Genscan
  chromosome:NCBIM37:11:27072113:27075899:1 // chr1 // 79
  // 73 // 19 // 24 // 0 /// GENSCAN00000020837 // ENSEMBL
  Prediction // cdna:Genscan
  chromosome:NCBIM37:4:8040558:8046315:1 // chr1 // 54 //
  73 // 13 // 24 // 0 /// GENSCAN00000035385 // ENSEMBL
  Prediction // cdna:Genscan
  chromosome:NCBIM37:3:4161125:4188920:-1 // chr1 // 12 //
  100 // 4 // 33 // 1" This should be 4 elements:
  fields[[1]] [1] "GENSCAN00000069186 // ENSEMBL Prediction
  // cdna:Genscan chromosome:NCBIM37:1:3018707:3044814:1 //
  chr1 // 100 // 100 // 33 // 33 // 0" [2]
  "GENSCAN00000048429 // ENSEMBL Prediction // cdna:Genscan
  chromosome:NCBIM37:11:27072113:27075899:1 // chr1 // 79
  // 73 // 19 // 24 // 0" [3] "GENSCAN00000020837 //
  ENSEMBL Prediction // cdna:Genscan
  chromosome:NCBIM37:4:8040558:8046315:1 // chr1 // 54 //
  73 // 13 // 24 // 0" [4] "GENSCAN00000035385 // ENSEMBL
  Prediction // cdna:Genscan
  chromosome:NCBIM37:3:4161125:4188920:-1 // chr1 // 12 //
  100 // 4 // 33 // 1" Each field has 9 parts, eg:
  fields[[1]] [[1]] [1] "GENSCAN00000069186" [2] "ENSEMBL
  Prediction" [3] "cdna:Genscan
  chromosome:NCBIM37:1:3018707:3044814:1" [4] "chr1" [5]
  "100" [6] "100" [7] "33" [8] "33" [9] "0" [[2]] [1]
  "GENSCAN00000048429" [2] "ENSEMBL Prediction" [3]
  "cdna:Genscan chromosome:NCBIM37:11:27072113:27075899:1"
  <snip>
}
\section{TODO}{
  - group columns into NCBI and ENSEMBL\cr - pick up
  miRNA's into a miRBase ID/description column
}
\author{
  Mark Cowley, 2008-07-21
}

