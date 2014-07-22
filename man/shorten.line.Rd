\name{shorten.line}
\alias{shorten.line}
\title{Shorten lines}
\usage{
shorten.line(coordinates, reduction)
}
\arguments{
  \item{coordinates}{length-4 vector of form c(x0, x1, y0,
  y1)}

  \item{reduction}{How much should length of the line be
  reduced by? Defaults to 0.05}
}
\value{
New coordinates of the line, with length reduced by amount
given by 'reduction'
}
\description{
Function to shorten arrow lengths, so that arrowheads are
always visible
}

