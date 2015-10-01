\name{mutual.information}
\alias{mutual.information}
\title{Calculate the association between two vectors using mutual information}
\usage{
mutual.information(dataset)
}
\arguments{
  \item{dataset}{Two column data.frame, with presence/absence of two species}
}
\value{
value of the mutual information between those species. Note this is necessarily less than the lower of their binary entropies.
}
\description{
Calculate the mututal information between two vectors
}

