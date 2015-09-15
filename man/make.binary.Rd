\name{make.binary}
\alias{make.binary}
\title{Convert a matrix or data.frame to binary}
\usage{
make.binary(dataset, threshold)
}
\arguments{
  \item{dataset}{a data.frame or matrix}

  \item{threshold}{the value above which values of
  'dataset' are set to one. Defaults to zero.}
}
\value{
Returns a numeric data.frame with same dimensions as
'dataset', but containing only zeros and ones.
}
\description{
Converts an abundance dataset to presence/absence.
}
\seealso{
cut
}

