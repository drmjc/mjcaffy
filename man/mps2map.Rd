\name{mps2map}
\alias{mps2map}
\title{Convert an Affymetrix MPS object into a 2column map, where each row contains
a key and a value; keys being the meta-probeset_id, and value being the
probeset_id's}
\usage{
  mps2map(x)
}
\arguments{
  \item{x}{an mps object. see import.mps}
}
\value{
  a map object: data.frame with 2 columns (metaprobeset_id
  and probeset_id)
}
\description{
  Convert an Affymetrix MPS object into a 2column map,
  where each row contains a key and a value; keys being the
  meta-probeset_id, and value being the probeset_id's
}
\examples{
\dontrun{
ex.mps3 <- mps2map(ex.mps)
head(ex.mps3)
}
}
\author{
  Mark Cowley, 2/1/2008
}

