\name{plot.NAA}
\alias{plot.NAA}
\title{Basic NAA plot}
\usage{
\method{plot}{NAA}(points, lines, draw.frequencies, add.key)
}
\arguments{
  \item{points}{data.frame giving point coordinates, of
  form given by or.points()}

  \item{lines}{data.frame giving line attributes, of form
  given by or.lines()}

  \item{draw.frequencies}{Should species occupancy rates be
  drawn as points? See 'notes' for possible values.}

  \item{add.key}{Should a key be added? See 'notes' for
  possible values.}
}
\value{
Draws a network diagram of species pairwise species
associations, based on strong +ve or -ve associations IDed
through odds ratios.
}
\description{
Basic NAA plot
}
\note{
Values of draw.frequencies are as follows: "none" = all
points the same size (default); "colour" = points are
coloured according to % occupied: darker = more common;
"size" = sized according to % occupied; "both" = higher %
occupied are larger and darker. Values of add.key are as
follows: "none" = no key; "objects" = points and line
widths; "species" = key to point labels; "both" = key to
species and point/line sizes
}

