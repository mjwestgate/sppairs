\name{or.lines}
\alias{or.lines}
\title{Line attributes}
\usage{
or.lines(dataset, threshold, reduce)
}
\arguments{
  \item{dataset}{output from pairwise.odds.ratios()}

  \item{thresold}{odds ratio used to determine groupings}

  \item{reduce}{if TRUE (default), looks for duplicated
  pairs of taxa, and returns the mean value for both pairs}
}
\description{
Function for setting the attributes of lines, to be drawn
between points given by or.points()
}

