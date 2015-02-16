\name{DRWPClass-GM-package}
\alias{DRWPClass-GM-package}
\alias{DRWPClass-GM}
\docType{package}
\title{
The DRW-GM method for disease classification
}
\description{
DRWPClass-GM performs pathway-based classifier construction and precise disease status prediction by joint analysis of genomic and metabolomic data and pathway topology.
}
\details{
\tabular{ll}{
Package: \tab DRWPClassGM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-02-11\cr
License: \tab GPL-2\cr
}
Very simple to use. \cr
Use \code{fit.DRWPClassGM} to train the classifier. \cr
Use \code{predict.DRWPClassGM} to predict the status of new samples.
}
\author{
Wei Liu, and Chunquan Li \cr
}
\references{
Liu, W., et al., Topologically inferring risk-active pathways toward precise cancer classification by directed random walk. 
  Bioinformatics, 2013. 29(17): p. 2169-77.
}
\keyword{ package }
\examples{
	data(GProf8511)
	data(GProf3325)
	data(MProf)
	data(pathSet)
	data(dGMGraph)
	fit <- fit.DRWPClassGM(xG=GProf8511$mRNA_matrix, yG.class1=GProf8511$normal, yG.class2=GProf8511$PCA,
	xM=MProf$Meta_matrix, yM.class1=MProf$normal, yM.class2=MProf$PCA, DEBUG=TRUE, pathSet=pathSet, 
	globalGraph=dGMGraph, testStatistic="t-test", classifier = "Logistic")
	predict.DRWPClassGM(object=fit, newx=GProf3325$mRNA_matrix[,c(GProf3325$normal,GProf3325$PCA)],	
	type = "class")
	evaluate.DRWPClassGM(object=fit, newx=GProf3325$mRNA_matrix, newy.class1=GProf3325$normal,
	newy.class2=GProf3325$PCA)
}
