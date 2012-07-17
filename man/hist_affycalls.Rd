\name{hist_affycalls}
\alias{hist_affycalls}
\title{Plot a histogram of affy P/M/A calls.}
\usage{
  hist_affycalls(calls, col = c("green", "orange", "red"),
    grid = TRUE, plot = TRUE, legend = TRUE,
    names.arg = colnames(calls), las = 2)
}
\arguments{
  \item{calls}{a data.frame, or expression set from
  mas5calls}

  \item{col}{the three colours corresponding to P, M, and
  A}

  \item{grid}{logical}

  \item{plot}{logical}

  \item{legend}{logical}

  \item{names.arg}{character vector}

  \item{las}{see ?par}
}
\value{
  none; generates a plot
}
\description{
  more correctly, this generates a stacked bar plot, where
  the counts of P, M and A are above each other, for each
  array
}
\examples{
\dontrun{
  par(mar=c(12,5,4,2), cex.lab=0.8)
  hist_affycalls(calls)
  grid()
  dev.off()
}
}
\author{
  Mark Cowley, 10 April 08
}

