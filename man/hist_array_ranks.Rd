\name{hist_array_ranks}
\alias{hist_array_ranks}
\title{histogram of array ranks}
\usage{
  hist_array_ranks(x, cols = NULL, abline = T, main = NULL,
    xlab = "rank (low to high)", auto.mfrow = T, ...)
}
\arguments{
  \item{x}{a matrix or data.frame containing the array
  data}

  \item{cols}{which columns to display? If NULL, use all
  the columns in x. If \code{x} is provided with N columns
  and only some are these are specified by \code{cols},
  then the ranks are still calculated on all the cols in
  \code{x}}

  \item{abline}{plot a horizontal line with the expected
  number of genes per rank?}

  \item{main}{either a single header, or 1 header per
  column. If NULL, then colnames will be used.}

  \item{xlab}{obvious}

  \item{auto.mfrow}{work out what par(mfrow=c(_,_)) should
  be used? see help for auto.mfrow.R}

  \item{\dots}{arguments passed to
  \code{\link[graphics]{hist}}}
}
\value{
  invisibly returns the ranks of x - same dimensions as x.
}
\description{
  plot each array as a histogram, to determine if the array
  has an unusual proportion of measurements that are
  consistently the highest (or lowest)
}
\details{
  It will plot upto 36 arrays in one panel. Any more than
  this would be too crowded. For each gene in 'x' (ie each
  row), determine the ranks across each array.
  Theoretically, all measurements on the array should be
  independent, so each array should have about the same
  number of genes with rank 1, or rank 2 or 3 etc... If an
  array has lots of genes that are rank 1, or rank
  \code{ncol(x)} then there is a problem with the array /
  normalisation
}
\author{
  Mark Cowley, 3 June 2005
}

