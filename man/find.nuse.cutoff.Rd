\name{find.nuse.cutoff}
\alias{find.nuse.cutoff}
\title{find a NUSE cutoff}
\usage{
  find.nuse.cutoff(x, file = NULL, nplots = 1,
    outlier.col = "red", type = c("plot", "outliers"))
}
\arguments{
  \item{x}{an rmaPLM object which has at least the standard
  errors stored}

  \item{file}{a pdf file to plot the nuse plots to (NULL
  means plot to dev.current())}

  \item{nplots}{If x has many samples, it may be useful to
  split the it into nplots}

  \item{outlier.col}{The colour to use for the outlying
  arrays}

  \item{type}{\dQuote{plot} means do the plots and
  invisibly return the outlier arrays; \dQuote{outliers}
  means don't do the plots and return the outlier arrays.}
}
\value{
  A \code{list} where each element is a \code{numeric
  vector} of indices into the columns/samples of x that are
  'outlier arrays'. The names of each element in the list
  indicate the threshold test that was applied.
}
\description{
  Take a \code{affyPLM} object obtained from fitting
  \code{\link[affyPLM]{fitPLM}}, and perform \code{NUSE}
  plots, identifying outlier based on various measure based
  on medians and IQR of the
  normalised-unscaled-standard-errors that NUSE calculates
}
\author{
  Mark Cowley, 8 Nov 2005
}

