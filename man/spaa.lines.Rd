\name{spaa.lines}
\alias{spaa.lines}
\title{Calculate line attributes}
\usage{
spaa.lines(dataset, threshold, reduce)
}
\arguments{
  \item{dataset}{output from spaa()}

  \item{thresold}{odds ratio used to determine groupings}

  \item{reduce}{if TRUE (default), looks for duplicated
  pairs of taxa, and returns the mean value for both pairs}
}
\description{
Function for setting the attributes of lines, to be drawn
between points given by spaa.points()
}

