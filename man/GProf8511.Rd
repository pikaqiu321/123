\name{GProf8511}
\alias{GProf8511}
\docType{data}
\title{
	Gene expression data from GSE8511.
}
\description{
	\code{GProf8511} contains gene expression profiles obtained from GSE8511 and the phenotypes of samples. 
}
\usage{data("GProf8511")}
\format{
  The format is:
List of 4
 $ mRNA_matrix: a 15728 x 41 matrix. Each row represents a gene, and each column represents a sample.
 $ normal     : An integer vector comprising the indexes of the Benign samples
 $ PCA        : An integer vector comprising the indexes of the PCA samples
 $ MET        : An integer vector comprising the indexes of the Mets samples
}
\examples{
data(GProf8511)
}
\keyword{datasets}
