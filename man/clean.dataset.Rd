\name{clean.dataset}
\alias{clean.dataset}
\title{Clean Dataset}
\usage{
clean.dataset(dataset, make.binary=TRUE, cutoff.min=0.1, cutoff.max=1, ...)
}
\arguments{
  \item{dataset}{a data.frame or matrix}

  \item{make.binary}{should the dataset be converted to binary? Defaults to TRUE}

  \item{cutoff.min}{minimum proportion of occupied rows for a species to be retained in the dataset.}

  \item{cutoff.max}{maximum proportion of occupied rows for a species to be retained in the dataset.}

  \item{...}{further information to pass to make.binary (i.e. to set a threshold)}
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

