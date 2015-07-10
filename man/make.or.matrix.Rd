\name{make.or.matrix}
\alias{make.or.matrix}
\title{Convert an asymmetric spaa result to a matrix}
\usage{
make.or.matrix(input)
}
\arguments{
  \item{input}{an object returned by spaa()}
}
\value{
a matrix containing the value of the asymmetric odds ratio between each pair of nodes
}
\description{
Converts the result from spaa to a square matrix. Useful for drawing results with library(circleplot) or similar.
}

