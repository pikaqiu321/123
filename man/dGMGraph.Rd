\name{dGMGraph}
\alias{dGMGraph}
\docType{data}
\title{
The global directed gene-metabolite pathway graph.
}
\description{
	The global directed gene-metabolite pathway graph constructed from the KEGG pathway dataset by the R package \code{iSubpathwayMiner}. 
}
\usage{data("dGMGraph")}
\format{
 An \code{igraph} R object.
}
\details{
There are 6386 nodes and 37939 edges in \code{dGMGraph}. Each node in the graph represents a gene or a metabolite.
The global pathway graph is used to evaluate the topological importances of genes by directed random walk.
}
\examples{
data(dGMGraph)
}
\keyword{datasets}
