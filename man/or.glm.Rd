\name{or.glm}
\alias{or.glm}
\title{Odds ratio calculation using regression}
\usage{
or.glm(dataset, complex)
}
\arguments{
  \item{dataset}{Two column data.frame containing species
  occurrence. Note that first column is species A and
  second species B, where A is the species subject to
  investigation (i.e. does the presence of B affect A?)}

  \item{complex}{logical, defaulting to FALSE. If true, or.glmer returns details of the model (slope and intercept) used to calculate each odds ratio, and the proportion of rows at which species B is present, as well as the odds ratio. Note that this may lead to incompatability with summary or plot methods as currently implemented - only recmomended for evaluation purposes.}

}
\value{
odds ratio between two species
}
\description{
Calculates asymmetric odds ratio using logistic regression.
}

