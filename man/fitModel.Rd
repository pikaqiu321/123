\name{fitModel}
\alias{fitModel}
\title{
Fitting model
}
\description{
This function trains classifier on the training set, and evaluates the classification performance of the classifier on the test set.
}
\usage{
fitModel(expr_training, expr_test, classifier)
}
\arguments{
  \item{expr_training}{
the training set.
}
  \item{expr_test}{
the test set
}
  \item{classifier}{
The method to train classifiers.
}
}
\value{
Fitted \code{"fitModel"} model object.
	\item{model}{Fitted model on the training set.}
	\item{eTestSet}{Classification performance of \code{model} on the test set.}
}
\seealso{
	\code{\link{fitModelGreedy}}
}
\examples{

}