\name{import.APT}
\alias{import.APT}
\title{Import a file produced by APT.}
\usage{
  import.APT(file, keep.first.column = FALSE,
    check.names = FALSE, ...)
}
\arguments{
  \item{file}{the file name}

  \item{keep.first.column}{logical: the first column can
  contain (IMO) useless row ID's.  if \code{TRUE}, then
  include them; if \code{FALSE}, exclude this first
  column.}

  \item{check.names}{logical: check the column names during
  the \code{\link{read.csv}}, or allow them to not strictly
  conform to \R's naming standards (\code{FALSE})}

  \item{\dots}{additional arguments passed to
  \code{\link{read.csv}}, or \code{\link{read.delim}}}
}
\value{
  a \code{data.frame} with row and col names, containing
  data from an APT analysis.
}
\description{
  Import a file produced by Affymetrix Power Tools (APT),
  (NOT debian's apt-get.org!)
}
\examples{
\dontrun{
rma <- import.APT("rma.summary.txt")
}
}
\author{
  Mark Cowley, 16 April 2008
}

