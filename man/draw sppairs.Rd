\name{draw.sppairs}
\alias{draw.sppairs}
\title{Draws basic sppairs plot}
\usage{
draw.sppairs(points, lines, draw.frequencies)
}
\arguments{
  \item{points}{data.frame giving point coordinates, of form given by or.points()}

  \item{lines}{data.frame giving line attributes, of form given by or.lines()}

  \item{draw.frequencies}{Should species occupancy rates be drawn as point sizes/colours? Defaults to TRUE}
}
\value{
Draws a network diagram of species pairwise species associations, based on strong +ve or -ve associations IDed through odds ratios. Typically called by plot.pairs().
}
\description{
Basic plot.
}

