\name{plot.spaa}
\alias{plot.spaa}
\title{Basic sppairs plot}
\usage{
plot.spaa(points, lines, draw.frequencies, add.key)
}
\arguments{
  \item{points}{data.frame giving point coordinates, of
  form given by or.points()}

  \item{lines}{data.frame giving line attributes, of form
  given by or.lines()}

  \item{draw.frequencies}{Should species occupancy rates be
  drawn as point sizes/colours? Defaults to TRUE}

  \item{add.key}{Should a key be added? See 'notes' for
  possible values.}
}
\value{
Draws a network diagram of species pairwise species
associations, based on strong +ve or -ve associations IDed
through SPAA.
}
\description{
Basic plot. Actually a wrapper function to allow plotting of a range of information from spaa()
}
\note{
Values of add.key are as follows: "none" = no key; "objects" = points and line widths; "species" = key to point labels; "both" = key to species and point/line sizes
}

