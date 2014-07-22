\name{pairwise.or.matrix}
\alias{pairwise.or.matrix}
\title{Pairwise odds ratio matrix}
\usage{
pairwise.or.matrix(or.object)
}
\arguments{
  \item{or.object}{data.frame containing columns listing
  the two species being compared, and the odds ratio
  between them}
}
\value{
A matrix
}
\description{
Internal function to convert the 'long' output from
pairwise.odds.ratios() into a matrix
}

