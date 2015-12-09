\name{spaa}
\alias{spaa}
\title{Species pairwise association analysis}
\usage{
spaa(x, method, asymmetric, ...)
}
\arguments{
  \item{x}{data.frame where each column is a species. Rows can be sites or visits; if visits, using or.glmer with a grouping factor might be a good idea.}

  \item{method}{Name of the function to be used to calculate the degree of association between two species.}

  \item{asymmetric}{If FALSE (the default), the pairwise association metric is assumed to be symmetric, and the function returns a data.frame with nrow(x)=choose(ncol(x), 2). If TRUE, the output has twice this many rows, with both directions of the association being calculated separately.}

  \item{...}{further information passed to 'method'}

}
\value{
A 3-column data.frame containing the names of the pair of species in question ($sp1, $sp2) the requested association metric ($value). If the specified function returns more than one value, then this function returns multiple columns (e.g. method="glmer.coef").
}
\description{
Calculates association metrics between all pairs of species in the dataset.
}
\note{
Be aware that this function can take a long time to run for large matrices. Methods for calculating pairwise associations in this package include:
  \describe{
    \item{fisher.test.pval}{Significance of association using Fisher's exact test}
    \item{glmer.coef}{Effect of species A on species B from mixed models; returns slope, se, z and p}
    \item{mutual.information}{Mutual information as defined in information theory}
    \item{or.asymmetric}{Asymmetric odds ratio using contingency tables}
    \item{or.glm}{Asymmetric odds ratio using logisitic regression}
    \item{or.glmer}{Asymmetric odds ratio using logisitic mixed models (using lme4)}
    \item{or.symmetric}{The default; symmetric odds ratio using contingency tables}
	}
}
\seealso{
\code{\link{plot.spaa}}
}

