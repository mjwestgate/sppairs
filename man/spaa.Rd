\name{spaa}
\alias{spaa}
\title{Pairwise odds ratio calculation}
\usage{
spaa(dataset, method, random.effect, rarity.cutoff)
}
\arguments{
  \item{dataset}{data.frame where each column is a species. Rows can be sites or visits; if visits, a grouping factor is a good idea.}

  \item{method}{value describing the method to be used. Should be one of "or.contingency", "glm", "glmer", "or.symmetric", "pearsons", "spearmans" or "mutual.info".}

  \item{random.effect}{Grouping variable for visits, if given. Passed to or.regression()}

  \item{rarity.cutoff}{Minimum proportion of occupied sites/visits for a species (column) to be included. Defaults to 0.1.}

}
\value{
A data.frame containing the requested association metric
}
\description{
Calculates association metrics between all pairs of (sufficiently common) species.
}
\seealso{
\code{\link{plot.spaa}}
}
