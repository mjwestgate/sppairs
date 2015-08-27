\name{spaa}
\alias{spaa}
\title{Pairwise odds ratio calculation}
\usage{
spaa(dataset, method, random.effect, rarity.cutoff)
}
\arguments{
  \item{dataset}{data.frame where each column is a species. Rows can be sites or visits; if visits, a grouping factor might be a good idea.}

  \item{method}{value describing the method to be used. Should be one of "or.contingency", "glm", "glmer", "or.symmetric", "pearson", "spearman" or "mutual.info".}

  \item{random.effect}{Grouping variable for visits, if given. Passed to or.regression()}

  \item{clean}{either a logical value (defaulting to TRUE) indicating whether to run clean.dataset(); or a list of commands to that function}

}
\value{
A data.frame containing the requested association metric
}
\description{
Calculates association metrics between all pairs of (sufficiently common) species.
}
\note{
Be aware that this function can take a long time to run for large matrices. Methods used to calculate pairwise associations are as follows:
  \describe{
    \item{or.contingency}{The default; calculates the asymmetric odds ratio using contingency tables}
    \item{glm}{Calculates the asymmetric odds ratio using logisitic regression}
    \item{glmer}{Calculates the asymmetric odds ratio using logisitic mixed models (using lme4)}
    \item{or.symmetric}{Calculates the symmetric odds ratio using contingency tables}
    \item{pearson}{Calculates Pearson's correlation coefficient}
    \item{spearman}{Calculates Spearman's correlation coefficient (rank correlation)}
    \item{mutual.info}{Calculates mutual information as defined in information theory}
	}
}
\seealso{
\code{\link{plot.spaa}}
}
