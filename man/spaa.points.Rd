\name{spaa.points}
\alias{spaa.points}
\title{Point locations}
\usage{
spaa.points(dataset, threshold, graph.function, ...)
}
\arguments{
  \item{dataset}{either an output from spaa(), or a list containing >1 of the same.}

  \item{thresold}{length-1 vector giving the odds ratio used to determine whether connections between pairs of taxa should be plotted. Both positive (x) and negative (1/x) associations are included in the output. Defaults to 3.}

  \item{graph.function}{if supplied, should be the name of a function (preferably from igraph) that can be used to determine the arrangement of points. Defaults to layout.fruchterman.reingold()}

  \item{...}{further commands to pass to the function given by graph.function (i.e. to determine point arrangement).}
}
\value{
A data.frame giving the frequency of occurrence and x, y coordinates for all included species. If many spaa() objects are supplied in a list, then the function returns a set of points that can be used to compare connections across all supplied objects.
}
\description{
Use igraph to calculate sensible point locations
}

