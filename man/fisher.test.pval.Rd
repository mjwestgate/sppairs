\name{fisher.test.pval}
\alias{fisher.test.pval}
\title{Calculate the significance of the association between two vectors using Fisher's exact test}
\usage{
fisher.test.pval(dataset)
}
\arguments{
  \item{dataset}{Two column data.frame, with presence/absence of two species}
}
\value{
Returns a value derived from the significance of the odds ratio between two species. Note that this version is experimental and has some unusual properties that might be subject to change. At present, the function runs two tests, to determine whether the odds ratio is significantly greater or less than one. It then returns 1-P, multiplied by -1 for negative associations. This means the resulting vector of pairwise assocations can be easily drawn (as it ranges from -1 to 1, with extreme values signifying strong -ve or +ve effects); but it is not a sensible statistic to report in this format.
}
\description{
Calculate the significance of the odds ratio between two vectors
}

