\name{spaa.points}
\alias{spaa.points}
\title{Point locations}
\usage{
spaa.points(dataset, threshold)
}
\arguments{
  \item{dataset}{output from spaa() (i.e. a
  list containing three components)}

  \item{thresold}{Odds ratio used to determine whether
  pairs of taxa should be plotted. Both positive (x) and
  negative (1/x) associations are included in the output.
  Defaults to 3.}
}
\value{
A data.frame giving the frequency of occurrence and x, y
coordinates for all included species.
}
\description{
Use igraph to calculate sensible point locations
}

