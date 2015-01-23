\name{plot.spaa}
\alias{plot.spaa}
\title{Network diagram of species pairwise associations}
\usage{
plot.spaa(object, plot.control, draw.frequencies)
}
\arguments{
  \item{object}{either an 'spaa' object (produced by command 'spaa'), or a list containing two data.frames giving giving point and line attributes (as returned by spaa.points and spaa.lines, in that order).}

  \item{plot.control}{allows detailed control of breaks, sizes and colours for points and lines. See details below.}

  \item{graph.function}{if supplied, should be the name of a function (preferably from igraph) that can be used to determine the arrangement of points. Passed to spaa.points; defaults to layout.fruchterman.reingold(). Note that supplying point and line values as a list will cause this function to be ignored.}

  \item{...}{further commands to pass to the function given by graph.function (i.e. to determine point arrangement).}
}
\value{
Draws a network diagram of species pairwise species associations, based on strong +ve or -ve associations identified using SPAA. This function also invisibly returns an object that can be passed to any of the functions that draw keys in sppairs (i.e. draw.spp.key, draw.line.key or draw.point.key).
}
\description{
Draws a network diagram of the results from species pairwise association analysis.
}
\note{
The plot options available have changed substantially from earlier versions of this package, and are now all manipulable using plot.control. If supplied, plot.control should be a named list containing any or all of the following components. 
  \describe{
    \item{threshold}{length-1 vector, giving the threshold for odds ratios to be included in the plot. Only lines >threshold or <(1/threshold) will be plotted.}
    \item{point.label}{vector. If length length-2, gives the colour for point edges and text (in that order); if length-1, both are set to the same value. Defaults to "grey30"}
    \item{point.breaks}{vector of breaks in species-level occupancy rates used to decide line colours Defaults to seq(0, 1, 0.2) if draw.frequencies=TRUE, else c(0, 1)}
    \item{point.cols}{vector containing colours, either named or in hexadecimal. Used with 'point.breaks' above; length should either equal 1 (same value for all points) or equal length(point.breaks)-1}
    \item{point.size}{a vector of sizes for points, with length equal to either one (all points the same size) or length(point.breaks)-1 (point size set in conjunction with point.breaks). Note that this is NOT passed to cex; point units are 1 percent of the size of the plot region. This tends to give larger points than those given by cex.}
    \item{line.breaks}{vector of breaks in species-level occupancy rates used to decide line colours Defaults to seq(0, 1, 0.2) if draw.frequencies=TRUE, else c(0, 1)}
    \item{line.cols}{vector containing colours, either named or in hexadecimal. Used with line.breaks above; length should equal length(line.breaks)-1}
    \item{line.widths}{if given, should provide line widths for each level defined by line.breaks, as for line.cols (above). Passed to lwd in lines().}
    \item{line.curvature}{length-1 vector, giving the proportion of the radius by which connecting lines should be diverted around intervening points. Defaults to 1.5; setting to zero will enforce straight lines between all vertices. This function is still in beta and may not lead to clean results in crowded plots.}
    \item{arrow.size}{length-1 numeric vector, used to set the arrowhead size, in the same units as point.size; defaults to 1.5}
    \item{arrow.angle}{length-1 numeric vector giving the angle that each edge of the arrowhead differs from the line direction, in degrees. Defaults to 20.}
    \item{text.size}{length-1 numeric vector passed to cex commands in all text options, allowing global scaling of label sizes.}

}
}
\seealso{
\code{\link{spaa}}
}
