\name{or.glmer}
\alias{or.glmer}
\title{Odds ratio calculation using mixed models}
\usage{
or.glmer(dataset, random.effect)
}
\arguments{
  \item{dataset}{Two column data.frame containing species
  occurrence. Note that first column is species A and
  second species B, where A is the species subject to
  investigation (i.e. does p/a of B affect A?)}

  \item{random.effect}{a factor to be used as random effect in glmer(lme4)}
}
\value{
value of the odds ratio between those species
}
\description{
Calculate asymmetric odds ratio using logistic regression.
Random effects must be specified, used e.g. to account for multiple
visits to the same site (using glmer in lme4). Also called by spaa().
}

