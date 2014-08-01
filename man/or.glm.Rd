\name{or.glm}
\alias{or.glm}
\title{Odds ratio calculation using regression}
\usage{
or.regression(dataset)
}
\arguments{
  \item{dataset}{Two column data.frame containing species
  occurrence. Note that first column is species A and
  second species B, where A is the species subject to
  investigation (i.e. does p/a of B affect A?)}
}
\value{
odds ratio between two species
}
\description{
Calculates asymmetric odds ratio using logistic regression.
}

