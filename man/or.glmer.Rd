\name{or.glmer}
\alias{or.glmer}
\title{Odds ratio calculation using mixed models}
\usage{
or.glmer(dataset, random.effect, complex)
}
\arguments{
  \item{dataset}{Two column data.frame containing species
  occurrence. Note that first column is species A and
  second species B, where A is the species subject to
  investigation (i.e. does p/a of B affect A?)}

  \item{random.effect}{a factor to be used as random effect in glmer(lme4)}

  \item{complex}{logical, defaulting to FALSE. If true, or.glmer returns details of the model (slope and intercept) used to calculate each odds ratio, and the proportion of rows at which species B is present, as well as the odds ratio. Note that this may lead to incompatability with summary or plot methods as currently implemented - only recommended for evaluation purposes.}
}
\value{
value of the odds ratio between those species
}
\description{
Calculate asymmetric odds ratio using logistic regression.
Random effects must be specified, used e.g. to account for multiple
visits to the same site (using glmer in lme4). Also called by spaa().
}

