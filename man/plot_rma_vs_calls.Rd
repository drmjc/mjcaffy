\name{plot_rma_vs_calls}
\alias{plot_rma_vs_calls}
\title{plot a distribution of RMA data, split by the affy P/M/A call.}
\usage{
  plot_rma_vs_calls(rma = NULL, calls = NULL,
    celfiles = NULL)
}
\arguments{
  \item{rma}{a \code{data.frame} of normalised data. If
  \code{NULL} then specify \code{celfiles}}

  \item{calls}{a \code{data.frame} of detection calls.  If
  \code{NULL} then specify \code{celfiles}}

  \item{celfiles}{a vector of cel files.}
}
\value{
  nonw. a plot is created
}
\description{
  either supply a vector of celfile names as full paths, or
  supply 2 data.frames (not eSet's) of RMA and MAS5 calls.
}
\author{
  Mark Cowley, 1/4/08
}

