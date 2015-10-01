\name{clean.dataset}
\alias{clean.dataset}
\title{Convert a dataset to binary, and remove unusually rare or common species}
\usage{
clean.dataset(x, make.binary=TRUE, cutoff.min=0.1, cutoff.max=1, ...)
}
\arguments{
  \item{x}{A data.frame or matrix}

  \item{make.binary}{Should the dataset be converted to binary? Defaults to TRUE}

  \item{cutoff.min}{Minimum proportion of occupied rows for a species to be retained in the dataset.}

  \item{cutoff.max}{Maximum proportion of occupied rows for a species to be retained in the dataset.}

  \item{...}{Further information to pass to make.binary (i.e. to set a threshold)}
}
\value{
Expansion to make.binary, that removes very rare and/or very common species from the dataset.
}
\description{
Converts a dataset to binary format and removes rare and/or common species.
}
\seealso{
cut
}

