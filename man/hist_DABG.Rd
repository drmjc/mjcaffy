\name{hist_DABG}
\alias{hist_DABG}
\title{Plot a histogram of DABG values}
\usage{
  hist_DABG(dabg,
    thresholds = c(1e-05, 1e-04, 0.001, 0.01, 0.1, 1),
    col = colour.step("green", "red", steps = length(thresholds)),
    grid = TRUE, plot = TRUE, legend = TRUE)
}
\arguments{
  \item{dabg}{the P-valuse from computing DABG calls.}

  \item{thresholds}{the P-value thresholds to use}

  \item{col}{the colours, 1 per threshold}

  \item{grid}{logical: overlay a grid?}

  \item{plot}{logical: create the plot?}

  \item{legend}{logical: add a legend to the plot,
  regarding the colour thresholds?}
}
\value{
  none.
}
\description{
  Plot a histogram of DABG values, thresholded at 10^-5,
  10^-4, 10^-3, 10^-2, 0.1 and 1.0
}
\author{
  Mark Cowley
}

