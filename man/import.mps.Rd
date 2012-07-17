\name{import.mps}
\alias{import.mps}
\title{Import an Affymetrix MPS file}
\usage{
  import.mps(file, ...)
}
\arguments{
  \item{file}{the file name}

  \item{\dots}{further arguments passed to read.delim (hint
  nlines=100)}
}
\value{
  a data.frame of 4 columns: probeset_id,
  transcript_cluster_id, probeset_list, probe_count
}
\description{
  Import an Affymetrix MPS file
}
\examples{
\dontrun{
ex.mps <- import.mps("~/data/Microarray Libraries/HuEx-1_0-st-v2/HuEx-1_0-st-v2.r2.dt1.hg18.core.mps", nrows=100)
head(ex.mps)
}
}
\author{
  Mark Cowley, 2/1/2008
}

