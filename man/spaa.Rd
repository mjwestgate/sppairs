\name{spaa}
\alias{spaa}
\title{Species pairwise association analysis}
\usage{
spaa(x, method, asymmetric, ...)
}
\arguments{
  \item{x}{data.frame where each column is a species. Rows can be sites or visits; if visits, using or.glmer with a grouping factor might be a good idea.}

  \item{method}{Name of the function to be used to calculate the degree of association between two species. Note that the specified function must return a single value.}

  \item{asymmetric}{If FALSE (the default), the pairwise association metric is assumed to be symmetric, and the function returns a data.frame with nrow(x)=choose(ncol(x), 2). If TRUE, the output has twice this many rows, with both directions of the association being calculated separately.}

  \item{...}{further information passed to 'method'}

}
\value{
A 3-column data.frame containing the names of the pair of species in question ($sp1, $sp2) the requested association metric ($value)
}
\description{
Calculates association metrics between all pairs of (sufficiently common) species.
}
\note{
Be aware that this function can take a long time to run for large matrices. Methods for calculating pairwise associations in this package include:
  \describe{
    \item{or.symmetric}{The default; Calculates the symmetric odds ratio using contingency tables}
    \item{or.asymmetric}{Calculates the asymmetric odds ratio using contingency tables}
    \item{or.glm}{Calculates the asymmetric odds ratio using logisitic regression}
    \item{or.glmer}{Calculates the asymmetric odds ratio using logisitic mixed models (using lme4)}
    \item{mutual.information}{Calculates mutual information as defined in information theory}
	}
}
\seealso{
\code{\link{plot.spaa}}
}

