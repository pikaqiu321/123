\name{getMetaWeight}
\alias{getMetaWeight}
\title{
	calculating the weights of metabolites.
}
\description{
	This function calculates the weights of metabolites in \code{metaProf} using Wilcoxon rank-sum test.
}
\usage{
getMetaWeight(metaProf, normMetaSmpl, diseaseMetaSmpl)
}
\arguments{
  \item{metaProf}{
a m x n matrix of metabolite expression measurements with m metabolites and n samples.
}
  \item{normMetaSmpl}{
a integer vector comprising the indexes of normal (class 1) samples in \code{metaProf}.
}
  \item{diseaseMetaSmpl}{
a integer vector comprising the indexes of disease (class 2) samples in \code{metaProf}.
}
}
\value{
The weighs of metabolites in \code{metaProf}.
}
\examples{
}