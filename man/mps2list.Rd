\name{mps2list}
\alias{mps2list}
\title{Convert an an Affymetrix MPS table into a named list, mapping affymtrix
meta-probeset_id's (the list element names) to the probeset_id's that make
up the meta-probeset (or the transcript cluster)}
\usage{
  mps2list(x)
}
\arguments{
  \item{x}{an mps object. see import.mps}
}
\value{
  a list of N entries, corresponding to N metaprobesets,
  where each element is a vector of probeset_id's.
}
\description{
  Convert an an Affymetrix MPS table into a named list,
  mapping affymtrix meta-probeset_id's (the list element
  names) to the probeset_id's that make up the
  meta-probeset (or the transcript cluster)
}
\examples{
\dontrun{
ex.mps2 <- mps2list(ex.mps)
}
}
\author{
  Mark Cowley, 2/1/2008
}

