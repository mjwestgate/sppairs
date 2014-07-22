\name{or.regression}
\alias{or.regression}
\title{Odds ratio calculation using regression}
\usage{
or.regression(dataset, random.effect)
}
\arguments{
  \item{dataset}{Two column data.frame containing species
  occurrence. Note that first column is species A and
  second species B, where A is the species subject to
  investigation (i.e. does p/a of B affect A?)}

  \item{random.effect}{# if specified, should be a factor
  to be used as random effect in glmer(lme4). If missing,
  glm() is used to calculate the odds ratio}
}
\value{
value of the odds ratio between those species
}
\description{
Calculate asymmetric odds ratio using logistic regression.
Random effects can be specified to account for multiple
visits to the same site (using glmer in lme4).
}

