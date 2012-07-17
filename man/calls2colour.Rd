\name{calls2colour}
\alias{calls2colour}
\title{Convert Affymetrix MAS5 calls into colours.}
\usage{
  calls2colour(x)
}
\arguments{
  \item{x}{If \code{x} is a \code{vector} of calls, then
  P/M/A = green/orange/red; else if \code{x} is a
  \code{matrix} of calls, then allP = green, noneP = red,
  else orange}
}
\value{
  a character vector or colours same length as
  \code{length(x)} or \code{nrow(x)}
}
\description{
  Convert Affymetrix MAS5 calls into colours.
}
\author{
  Mark Cowley, 2009-05-15
}

