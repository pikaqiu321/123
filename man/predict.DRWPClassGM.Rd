\name{predict.DRWPClassGM}
\alias{predict.DRWPClassGM}
\title{
Make predictions from a "DRWPClassGM" object.
}
\description{
This functions predicts class labels or probabilities of unknown samples from a fitted \code{"DRWPClassGM"} object.
}
\usage{
predict.DRWPClassGM(object, newx, type = c("class", "probability"))
}
\arguments{
  \item{object}{
Fitted \code{"DRWPClassGM"} model object.
}
  \item{newx}{
A matrix with variables to predict.
}
  \item{type}{ \code{"class"} gives the predicted classes; \code{"probability"} gives the predictions.}
}
\value{
the predictions.
}
\seealso{
\code{\link{fit.DRWPClassGM}}
}
\examples{
	data(GProf8511)
	data(GProf3325)
	data(MProf)
	data(pathSet)
	data(dGMGraph)
	fit <- fit.DRWPClassGM(xG=GProf8511$mRNA_matrix, yG.class1=GProf8511$normal, yG.class2=GProf8511$PCA,
                           xM=MProf$Meta_matrix, yM.class1=MProf$normal, yM.class2=MProf$PCA, DEBUG=TRUE, 
                           pathSet=pathSet,	globalGraph=dGMGraph, testStatistic="t-test", classifier = "Logistic", 
                           normalize = TRUE, nFolds = 5, numTops=50, iter = 1, Gamma=0.7, Alpha = 0.5)
	predict.DRWPClassGM(object=fit, newx=GProf3325$mRNA_matrix[,c(GProf3325$normal,GProf3325$PCA)],	
                        type = "class")
}