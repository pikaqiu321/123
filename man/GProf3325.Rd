\name{GProf3325}
\alias{GProf3325}
\docType{data}
\title{
	Gene expression data from GSE3325.
}
\description{
	\code{GProf3325} contains gene expression profiles obtained from GSE3325 and the phenotypes of samples. 
}
\usage{data("GProf3325")}
\format{
  The format is:
List of 4
 $ mRNA_matrix: a 15728 x 19 matrix. Each row represents a gene, and each column represents a sample. 
 $ normal     : An integer vector comprising the indexes of the Benign samples
 $ PCA        : An integer vector comprising the indexes of the PCA samples
 $ MET        : An integer vector comprising the indexes of the Mets samples
}
\examples{
data(GProf3325)
}
\keyword{datasets}
