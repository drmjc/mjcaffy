\name{readCEL}
\alias{readCEL}
\title{Function to import an Affymetrix CEL file}
\usage{
  readCEL(x, type = NULL, skip = 24, numlines = NULL)
}
\arguments{
  \item{x}{The file name. Can be a vector of filenames.}

  \item{type}{MG_U74Av2 | HG_Focus | RAE230A | NULL}

  \item{skip}{The number of lines to skip (if known)}

  \item{numlines}{The number of lines in the [INTENSITY]
  section to read (if known)}
}
\value{
  if x is only 1 file, then a \code{data.frame} is returned
  with 5 columns: "X", "Y", "mean", "sd", "npixels" if
  length(x) > 1 then a \code{list} of \code{data.frames} is
  returned.
}
\description{
  The file can be unzipped, or gzipped. Since we are only
  interested in the [INTENSITY] section, we need to find
  where to start and stop reading the file. If the type is
  one of the known array types, please specify this in the
  'type' parameter. If left as NULL, then each file will be
  searched for when to start and stop reading which takes
  longer.
}
\author{
  Mark Cowley, 28 Oct 2005
}

