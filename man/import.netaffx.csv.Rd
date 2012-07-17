\name{import.netaffx.csv}
\alias{import.netaffx.csv}
\title{import Affymetrix NetAffx csv file}
\usage{
  import.netaffx.csv(file, ...)
}
\arguments{
  \item{file}{the path to a NetAffx csv file}

  \item{\dots}{other arguments passed to
  \code{\link{read.csv}}}
}
\value{
  a \code{data.frame}
}
\description{
  import a file from NetAffx, designed to work with
  annotation files that have comment rows, and are csv's.
  Warning, this can take a long time! 41s to read the HuEx
  na32 transcript.csv from an SSD on a 2010 MBP, or 192s to
  read the HuEx na32 probeset.csv file
}
\author{
  Mark Cowley, 2008-05-27
}

