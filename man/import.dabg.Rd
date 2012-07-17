\name{import.dabg}
\alias{import.dabg}
\title{Import DABG results}
\usage{
  import.dabg(file = "dabg.summary.txt", ...)
}
\arguments{
  \item{file}{the path to the DABG file}

  \item{\dots}{arguments passed to
  \code{\link{import.APT}}}
}
\value{
  a \code{data.frame} containing DABG P-values.
}
\description{
  Import a \dQuote{dabg.summary.txt} file from running the
  Affymetrix Detected Above Background method (DABG),
  skipping the header
}
\author{
  Mark Cowley, 11/3/08
}
\seealso{
  \code{\link{dabg2calls}} \code{\link{import.APT}}
}

