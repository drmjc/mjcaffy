\name{parse.affy.ST.array.annotation}
\alias{parse.affy.ST.array.annotation}
\title{parse an Affymetrix Gene or Exon ST transcript.csv file}
\usage{
  parse.affy.ST.array.annotation(csv.file, mps.file = NULL,
    probes = NULL, char.limit = 512, create.files = TRUE,
    verbose = TRUE)
}
\arguments{
  \item{csv.file}{path to the transcript.csv file}

  \item{mps.file}{Useful Exon arrays to subset the
  transcript-level probesets to the core, extended or full.
  The mps file specifies which probesets are in each
  category.  This should be ignored for Gene Arrays.}

  \item{probes}{a character vector of probe ID's. If
  \code{NULL} (the default), then all probes in the csv
  file will be returned. if not \code{NULL}, then only
  those probes will be returned, in the same order as
  probes.}

  \item{char.limit}{restrict the length of the text in each
  column to something manageble}

  \item{create.files}{logical: create the tsv and RDa
  files?}

  \item{verbose}{logical: verbose output?}
}
\value{
  For a given an input file, example:
  \code{csv.file="HuGene-1_0-st-v1.na26.hg18.transcript.csv"}\cr
  These files are created:\cr
  \dQuote{HuGene-1_0-st-v1.na26.hg18.transcript.Rda.gz},
  and\cr
  \dQuote{HuGene-1_0-st-v1.na26.hg18.gene.info.Rda.gz},
  and\cr \dQuote{HuGene-1_0-st-v1.na26.hg18.gene.info.tsv},
  and\cr \code{gene.info} is invisibly returned
}
\description{
  parse an Affymetrix Gene or Exon ST transcript.csv file,
  and:\cr - create gene.info.raw (essentially an R version
  of the full csv file)\cr - create gene.info gene.info
  contains 6 columns:\cr ProbeSetID\cr geneSymbol (the 2nd
  "//" seperated term in gene_assignment)\cr Description
  (the 3nd "//" seperated term in gene_assignment)\cr
  GenBank ID (the 1st "//" seperated term in
  mrna_assignment, \code{NA}ing out non-GenBanks (ie
  ENS...))\cr gene_assignment\cr mrna_assignment\cr
}
\details{
  Specifying which probes:\cr Sometimes there are missing
  probes that are in the PGF/CLF file, but not in the CSV
  file. To be sure, provide a vector of probe ID's for the
  expected probeset ID's. For those probes that are missing
  annotation, default annotations will be reported.
}
\author{
  Mark Cowley, 2008-07-24
}

