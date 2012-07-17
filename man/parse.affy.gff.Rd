\name{parse.affy.gff}
\alias{parse.affy.gff}
\title{parse.affy.gff}
\usage{
  parse.affy.gff(gff, key)
}
\arguments{
  \item{gff}{A gff object with at least a \code{$group}
  attribute}

  \item{key}{Unknown}
}
\value{
  Undocumented return value
}
\description{
  Affy GFF files are "nearly" GTF files in that the
  9th/last column contains key/value pairs. nearly, because
  GTF files must begin with the two keys of "gene_id" and
  "transcript_id". The pairs are "; " seperated, last field
  may not have the ";"
}
\details{
  The 9th field in a gff table contains lots of semicolon
  seperated key/value pairs. If you'd like to retrieve the
  value, associated with a particular key, then this method
  is the goods. eg: "number_independent_probes 3;
  probeset_id 4027908; exon_cluster_id 1056646;"
}
\author{
  Mark Cowley, 2/1/08
}

