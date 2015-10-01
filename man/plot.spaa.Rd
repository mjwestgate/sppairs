\name{plot.spaa}
\alias{plot.spaa}
\title{Network diagram of species pairwise associations}
\usage{
plot.spaa(object, plot.control, draw.frequencies)
}
\arguments{
  \item{object}{a data.frame, in which the first two columns give the names of each pair of nodes, and the third gives the degree of association between those node. This is the format of data returned by spaa().}

  \item{...}{further commands passed to plot.igraph}
}
\value{
Draws a network diagram of species pairwise species associations, based on strong +ve or -ve associations identified using the selected association function.
}
\description{
A wrapper function to call igraph in the correct manner to plot spaa() results.
}

\seealso{
\code{\link{spaa}}
}
