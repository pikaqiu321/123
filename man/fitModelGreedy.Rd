\name{fitModelGreedy}
\alias{fitModelGreedy}
\title{
Fitting model using a greedy search
}
\description{
This function selects pathway features from the top \code{numTops} pathways using a greedy search and returns the Fitted \code{"fitModel"} model object with best AUC.
}
\usage{
fitModelGreedy(expr_training, expr_test, TValue, classifier = "Logistic", numTops = 50)
}
\arguments{
  \item{expr_training}{
	The training set.
}
  \item{expr_test}{
	The test set.
}
  \item{TValue}{
	The sorted t statistic of pathway features.
}
  \item{classifier}{
	The method to train classifiers.
}
  \item{numTops}{
	The number of pathway features used for feature selection. 
}
}
\value{
	\item{model}{Fitted \code{"fitModel"} model object.}
	\item{AUC}{The performance (AUC) of the classifier on \bold{feature selection set}.}
	\item{Accuracy}{The performance (Accuracy) of the classifier on \bold{feature selection set}.}
	\item{pathFeatures}{The selected pathway features to build the classifier.}
}
\seealso{
 \code{\link{fitModel}}
}
\examples{
}