\name{glmer.coef}
\alias{glmer.coef}
\title{Calculate pairwise association between species using mixed models}
\usage{
glmer.coef(dataset, random.effect)
}
\arguments{
  \item{dataset}{Two column data.frame containing species
  occurrence. Note that first column is species A and
  second species B, where A is the species subject to
  investigation (i.e. does p/a of A affect B?)}

  \item{random.effect}{a factor to be used as random effect in glmer(lme4)}

}
\value{
Returns information on the effect of species B on A using summary(model)$coefficients[2, ]. Gives the slope, standard error, z score and p value. 
}
\description{
Calculate the effect of one species on a second species using logistic regression.
Random effects must be specified, and can be used to account for multiple
visits to the same site (for example). Calculations of mixed models are from glmer in lme4. 
}

