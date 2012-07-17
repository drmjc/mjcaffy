\name{cdfName}
\alias{cdfName}
\alias{cdfName,character-method}
\title{Extract the CDF file name from a CEL file header.}
\arguments{
  \item{object}{either an \code{\link[affy]{AffyBatch}}, or
  character vector of filenames}
}
\value{
  a vector of cdf names for each CEL file.
}
\description{
  affy provides \code{cdfName,AffyBatch}, this provides
  \code{cdfName,character}, where a character vector of
  file paths can be used to obtain the CDF name embedded in
  the CEL file header.
}
\examples{
\dontrun{
files <- dir("/Volumes/Volumes/PWBC/private/projects/ShaneGrey/CEL file repsitory/CEL/", full.names=TRUE, pattern="AW.*CEL")[1:3]
cdfName(files)
}
}
\author{
  Mark Cowley, 2008-05-26
}

