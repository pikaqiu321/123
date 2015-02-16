\name{evaluate.DRWPClassGM}
\alias{evaluate.DRWPClassGM}
\title{
Evaluate predictions from a "DRWPClassGM" object.
}
\description{
This functions evaluates the classification performance of a fitted \code{"DRWPClassGM"} object. 
}
\usage{
evaluate.DRWPClassGM(object, newx, newy.class1, newy.class2)
}
\arguments{
  \item{object}{
Fitted \code{"DRWPClassGM"} model object.
}
  \item{newx}{
A matrix with variables to predict.
}
  \item{newy.class1}{
a integer vector comprising the indexes of class 1 samples in \code{newx}.
}
  \item{newy.class2}{
a integer vector comprising the indexes of class 2 samples in \code{newx}.
}
}
\value{The classification performance of the predictions.}
\seealso{
\code{\link{predict.DRWPClassGM}}
}
\examples{
	data(GProf8511)
	data(GProf3325)
	data(MProf)
	data(pathSet)
	data(dGMGraph)
	fit <- fit.DRWPClassGM(xG=GProf8511$mRNA_matrix, yG.class1=GProf8511$normal, 
						   yG.class2=GProf8511$PCA, xM=MProf$Meta_matrix, yM.class1=MProf$normal, 
						   yM.class2=MProf$PCA,	DEBUG=TRUE, pathSet=pathSet, globalGraph=dGMGraph, 
						   testStatistic="t-test", classifier = "Logistic", normalize = TRUE, 
						   nFolds = 5, numTops=50, iter = 1, Gamma=0.7, Alpha = 0.5)
	evaluate.DRWPClassGM(object=fit, newx=GProf3325$mRNA_matrix, newy.class1=GProf3325$normal, 
						 newy.class2=GProf3325$PCA)
}