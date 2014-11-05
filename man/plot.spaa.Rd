\name{plot.spaa}
\alias{plot.spaa}
\title{Basic sppairs plot}
\usage{
plot.spaa(object, plot.control, draw.frequencies)
}
\arguments{
  \item{object}{either an 'spaa' object (produced by command 'spaa'), or a (preferably named) list containing two data.frames giving giving point and line attributes.}

  \item{plot.control}{allows detailed control of breaks, sizes and colours for points and lines. See details below.}

  \item{draw.frequencies}{Should species occupancy rates be
  drawn as point sizes/colours? Defaults to TRUE. Note that this is ignored if new point attributes are supplied using plot.control}

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
The plot options available have changed substantially from earlier versions of this package, and are now all manipulable using plot.control. If supplied, plot.control should be a named list containing any or all of the following components. 
  \describe{
    \item{point.label}{length-1 character vector giving the colour for point edges and text. Defaults to "grey30"}
    \item{point.breaks}{vector of breaks in species-level occupancy rates used to decide line colours Defaults to seq(0, 1, 0.2) if draw.frequencies=TRUE, else c(0, 1)}
    \item{point.cols}{vector containing colours, either named or in hexadecimal. Used with 'point.breaks' above; length should equal length(point.breaks)-1}
    \item{point.size}{either a single value giving the size of all points, or a range giving the min and max values, which are then scaled according to occupancy rates.}
    \item{line.breaks}{vector of breaks in species-level occupancy rates used to decide line colours Defaults to seq(0, 1, 0.2) if draw.frequencies=TRUE, else c(0, 1)}
    \item{line.cols}{vector containing colours, either named or in hexadecimal. Used with line.breaks above; length should equal length(line.breaks)-1}
    \item{line.widths}{if given, should provide line widths for each level defined by line.breaks, as for line.cols (above).}
    \item{text.size}{length-1 numeric vector passed to cex commands in all text options, allowing global scaling of label sizes.}
    \item{key.placement}{If given, should be a matrix with four columns giving the left, right, bottom and top proportions of the plot window to be used for each section of the plot. Further, the rows (ip to four) should be named with the kind of plot in each section; namely one of 'network' (plotted alone by default) and keys to "species", "points" and/or "lines". This matrix is passed directly to split.screen; help for that command may assist with understanding use of this matrix.}
}
}