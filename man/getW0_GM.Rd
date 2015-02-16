\name{getW0_GM}
\alias{getW0_GM}
\title{
calculating the weights of genes and metabolites
}
\description{
This function calculates the initial weights of genes and metabolites in \code{globalGraph}. The initial weights are used to start the directed random walk on \code{globalGraph}.
}
\usage{
getW0_GM(geneWeight, metaWeight = NULL, globalGraph)
}
\arguments{
  \item{geneWeight}{
The weighs of genes.
}
  \item{metaWeight}{
The weighs of metabolites.
}
  \item{globalGraph}{
An \code{igraph} R object contains the global directed gene-metabolite pathway graph.
}
}
\value{
The initial weights of genes and metabolites in \code{globalGraph}
}
\seealso{
\code{\link{getMetaWeight}}, \code{\link{calSAMScore}}, \code{\link{caltScore}}
}
\examples{

}