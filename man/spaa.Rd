\name{spaa}
\alias{spaa}
\title{Pairwise odds ratio calculation}
\usage{
spaa(dataset, method, random.effect, rarity.cutoff)
}
\arguments{
  \item{dataset}{data.frame where each column is a species. Rows can be sites or visits; if visits, a grouping factor is a good idea.}

  \item{method}{value describing the method to be used. Should be one of "or.contingency", "glm", "glmer", "or.symmetric", "pearsons" or "spearmans".}

  \item{random.effect}{Grouping variable for visits, if given. Passed to or.regression()}

  \item{rarity.cutoff}{Minimum proportion of occupied sites/visits for a species (column) to be included. Defaults to 0.1.}

}
\value{
A list containing the components: [[1]] the frequency at which species occur at sites; [[2]] pairwise odds ratios in long format; [[3]] pairwise odds ratios in wide format
}
\description{
Calculates a number of odds ratios for all (sufficiently common) species. Basically a wrapper function for
or.contingency() and or.regression().
}
\seealso{
\code{\link{plot.spaa}}
}
