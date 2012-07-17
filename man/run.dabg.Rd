\name{run.dabg}
\alias{run.dabg}
\title{run Affymetrix DABG algorithm}
\usage{
  run.dabg(cel.files = NULL, lib.dir = NULL,
    interactive = FALSE)
}
\arguments{
  \item{cel.files}{a character vector of cel file names}

  \item{lib.dir}{the location of the PGF/CLF files; or set
  \code{NULL} & specify \code{interactive=TRUE}}

  \item{interactive}{logical: if \code{TRUE}, the ask user
  for the CEL file dir, PGF, CLF, BGP files}
}
\value{
  the DABG calls as a \code{data.frame}
}
\description{
  Perform DABG method to determine P-values for each
  ProbeSet x array
}
\section{TODO}{
  use \code{\link{Sys.which}}
}
\author{
  Mark Cowley, 29/1/08
}

