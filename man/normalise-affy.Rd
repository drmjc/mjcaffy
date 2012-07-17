\name{norm.dchip}
\alias{norm.dchip}
\alias{norm.mas5}
\alias{norm.mas5.slow}
\alias{norm.rma}
\title{normalize Affymetrix arrays}
\usage{
  norm.dchip(cel.data, ...)

  norm.rma(cel.data, ...)

  norm.mas5.slow(cel.data, scale = 500, ...)

  norm.mas5(cel.data, scale = 500, ...)
}
\arguments{
  \item{cel.data}{an \code{AffyBatch} object}

  \item{\dots}{arguments passed to sub-calls}

  \item{scale}{Value at which all arrays will be scaled
  to.}
}
\value{
  An object of class \code{ExpressionSet}
}
\description{
  normalize Affymetrix arrays via dChip, rma, mas5
}
\author{
  Mark Cowley, 27 April 2005
}

