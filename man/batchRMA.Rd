\name{batchRMA}
\alias{batchRMA}
\title{RMA normalize in batches}
\usage{
  batchRMA(files, batch)
}
\arguments{
  \item{files}{a vector of CEL filenames}

  \item{batch}{a vector that identifies the batches. see
  Details.}
}
\value{
  a \code{data.frame} of RMA normalised data;
}
\description{
  Normalize a set of arrays, using RMA in batches according
  to the date that they were processed.
}
\details{
  Each array's batch can be specified via anything that
  uniquely identifies arrays from the same batch, eg
  numerical ID, a character name, or a character date
}
\note{
  further normalisation may be (usually is) required
  between batches, eg quantile norm.
}
\examples{
\dontrun{
rma.raw <- batchRMA(files, c(1,1,1,2,2,2,3,3,3))
rma <- normalizeQuantiles(rma.raw)
}
}
\author{
  Mark Cowley, 10 April 08
}

