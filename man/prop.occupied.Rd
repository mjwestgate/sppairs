\name{prop.occupied}
\alias{prop.occupied}
\title{Proportion of occupied rows}
\usage{
prop.occupied(dataset, random.effect)
}
\arguments{
  \item{dataset}{a data frame or matrix containing
  presence/absence site x species data}

  \item{random.effect}{if given, function calculates the
  proportion from xtab of levels of this variable.}
}
\value{
A vector with length = ncol(dataset), giving the proportion
of sites/visits occupied by each species.
}
\description{
Calculate the proportion of rows occupied (where rows
represent either visits or sites)
}

