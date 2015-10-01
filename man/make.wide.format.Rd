\name{make.wide.format}
\alias{make.wide.format}
\title{Convert an asymmetric spaa result to a matrix}
\usage{
make.wide.format(input)
}
\arguments{
  \item{input}{a data.frame, in which the first two columns give the names of each pair of nodes, and the third gives the degree of association between those node. This is the format of data returned by spaa().}
}
\value{
a matrix containing the value of the association function between each pair of nodes
}
\description{
Converts the result from spaa to a square matrix. Useful for drawing results with library(circleplot) or similar.
}

