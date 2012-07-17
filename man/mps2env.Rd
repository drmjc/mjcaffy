\name{mps2env}
\alias{mps2env}
\title{Convert an Affymetrix MPS object into a hash table (an R environment); The
keys being the meta-probeset_id, and value being the probeset_id's}
\usage{
  mps2env(x)
}
\arguments{
  \item{x}{an mps object. see import.mps}
}
\value{
  a hash table/R environment
}
\description{
  Convert an Affymetrix MPS object into a hash table (an R
  environment); The keys being the meta-probeset_id, and
  value being the probeset_id's
}
\author{
  Mark Cowley, 2/1/2008
}

