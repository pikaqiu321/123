\name{MProf}
\alias{MProf}
\docType{data}
\title{
	Metabolomic data
}
\description{
	The metabolomic profiles, including 16 Benign, 12 PCA and 14 Mets, were obtained from studies of Sreekumar et al.
}
\usage{data("MProf")}
\format{
  The format is:
List of 4
 $ Meta_matrix: a 145 x 42 matrix. Each row represents a metabolite, and each column represents a sample. 
 $ normal     : An integer vector comprising the indexes of the Benign samples
 $ PCA        : An integer vector comprising the indexes of the PCA samples
 $ MET        : An integer vector comprising the indexes of the Mets samples
}
\examples{
data(MProf)
}
\keyword{datasets}
