\name{getPathActivity}
\alias{getPathActivity}
\title{
inferring pathway activity
}
\description{
This function infers pathway activities using the DRW-GM method.
}
\usage{
getPathActivity(x, pathSet, w, vertexZP)
}
\arguments{
  \item{x}{
a p x n matrix of gene expression measurements with p genes and n samples.
}
  \item{pathSet}{
	A list of pathways. Each pathway is represented as a vector of pathway member genes and metabolites.
}
  \item{w}{
	A numerical vector containing the topological weights of nodes in \code{globalGraph}.
}
  \item{vertexZP}{
	A p x 2 matrix which contains the t statistic and p-value (d statistic and q-value) of p genes in \code{x}.
}
}
\details{
For each pathway, we selected the differential genes in the pathway to infer pathway activity. The expression values of genes are weighted by their topological weights obtained from directed random walk on the global directed gene-metabolite graph. Pathway activity a(Pj) is a synthetic expression value vector of Pj across the samples. Pathway activity inference transforms the gene expression profiles into pathway activity profiles, which are then used to train a classifier.
}
\value{
	\item{pathwayActivity}{The pathway activities of pathways in \code{pathSet}}
	\item{sigGenes}{The genes used to infer the pathway activity.}
}
\references{
Liu, W., et al., Topologically inferring risk-active pathways toward precise cancer classification by directed random walk. 
  Bioinformatics, 2013. 29(17): p. 2169-77.
}
\author{
Wei Liu
}
\examples{
	data(dGMGraph)
	data(pathSet)
	data(GProf8511)
	tScore <- caltScore(GProf8511$mRNA_matrix, GProf8511$normal, GProf8511$PCA)
	
	vertexs <- V(dGMGraph)
	p0 <- runif(length(vertexs), min = 0, max = 1)
	names(p0) <- vertexs$name
	p0 <- p0/sum(p0)
	vertexWeight <- DRW(igraphM = dGMGraph, p0=p0, EdgeWeight=FALSE, gamma = 0.3)
	names(vertexWeight) <- names(p0)

	pA <- getPathActivity(GProf8511$mRNA_matrix[ ,c(GProf8511$normal, GProf8511$PCA)],
						  pathSet, vertexWeight, tScore)
}