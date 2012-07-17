\name{dabg2calls}
\alias{dabg2calls}
\title{Convert DABG P-values into P/M/A calls}
\usage{
  dabg2calls(dabg, Pthresh = 1e-05, Mthresh = 0.001)
}
\arguments{
  \item{dabg}{a data.frame of dabg calls. see ?import.APT}

  \item{Pthresh}{the p-values at which to threshold the
  DABG p-values}

  \item{Mthresh}{the p-values at which to threshold the
  DABG p-values}
}
\value{
  a data.frame, with same dim, and dimname as dabg, but
  with values of "P", "M", "A"
}
\description{
  Convert DABG P-values into P/M/A calls, by
  hard-thresholding at two threshold, one for "Present",
  and "Marginal"
}
\author{
  Mark Cowley, 2008-05-26
}

