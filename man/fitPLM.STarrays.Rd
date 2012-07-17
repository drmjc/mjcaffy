\name{fitPLM.STarrays}
\alias{fitPLM.STarrays}
\title{fitPLM model on Affymetrix ST arrays}
\usage{
  fitPLM.STarrays(object, ..., subset = NULL,
    max.probes = 1000, verbose = TRUE)
}
\arguments{
  \item{object}{see fitPLM. My testing assumes that object
  is an AffyBatch.}

  \item{\dots}{see fitPLM. My testing assumes that object
  is an AffyBatch.}

  \item{subset}{see fitPLM. My testing assumes that object
  is an AffyBatch.  subset is over-written by this
  function.}

  \item{max.probes}{determine the largest number of probes
  allowed in a probeset.}

  \item{verbose}{TRUE/FALSE to print a message about how
  many probes are removed.}
}
\value{
  the outpur from running fitPLM on a subset of probesets
}
\description{
  Affymetrix ST arrays have some probesets with thousands
  of probes which makes fitPLM REALLY slow, and prone to
  running out of memory. Limit the probesets to those with
  at most max.probes. Most of these probes that are
  excluded are control probes, but not all.
}
\details{
  Using a 4x2.66 GHz Mac Pro with 5 GB RAM, and 6 MoGene
  arrays, I find max.probes [200,1000] gives similar run
  time. limiting to 100 probes = 74 s 200 probes = 69 s 500
  probes = 67 s 1000 probes = 80 s all probes = > 4min
  before I get a memory allocation error. My 1.8 GHz
  MacBook Pro can run probesets with < 1000 probes in 134s
  which is acceptable.
}
\author{
  Mark Cowley, 2009-06-24
}

