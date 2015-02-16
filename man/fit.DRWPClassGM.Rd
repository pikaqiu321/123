\name{fit.DRWPClassGM}
\alias{fit.DRWPClassGM}
\title{
Train a classifier using DRW-GM
}
\description{
DRW-GM is a disease classification method which performs pathway-based classifier construction and precise disease status prediction by joint analysis of genomic and metabolomic data and pathway topology.
}
\usage{
fit.DRWPClassGM(xG, yG.class1, yG.class2, xM = NULL, yM.class1 = NULL, 
                yM.class2 = NULL, DEBUG = FALSE, pathSet, globalGraph,
                testStatistic = c("t-test", "SAM"), classifier = "Logistic", 
                normalize = TRUE, nFolds = 5, numTops = 50, iter = 1, 
                Gamma = 0.7, Alpha = 0.5, fdr.output = 0.2)
}
\arguments{
  \item{xG}{
	a p x n matrix of gene expression measurements with p genes and n samples.
}
  \item{yG.class1}{
	a integer vector comprising the indexes of class 1 samples in \code{xG}.
}
  \item{yG.class2}{
	a integer vector comprising the indexes of class 2 samples in \code{xG}.
}
  \item{xM}{
	a m x n matrix of metabolite expression measurements with m metabolites and n samples.
}
  \item{yM.class1}{
	a integer vector comprising the indexes of class 1 samples in \code{xM}.
}
  \item{yM.class2}{
	a integer vector comprising the indexes of class 2 samples in \code{xM}.
}
  \item{DEBUG}{
	Logical. Should debugging information be plotted.
}
  \item{pathSet}{
	A list of pathways. Each pathway is represented as a vector of pathway member genes and metabolites.
}
  \item{globalGraph}{
	 An \code{igraph} R object contains the global directed gene-metabolite pathway graph.
}
  \item{testStatistic}{
	The test method used to identify differential genes. For \code{testStatistic="t-test"}, function \code{caltScore} is used. For \code{testStatistic="SAM"}, function \code{calSAMScore} is used.
}
  \item{classifier}{
	The method to train classifiers. The default is "Logistic". To use other methods, such as "libsvm", one should install the corresponding package in Weka.
}
  \item{normalize}{
	Logical flag for \code{xG} and \code{xM} standardization, prior to fitting the model. Default is \code{normalize=TRUE}.
}
  \item{nFolds}{
	The number of folds to split \code{xG}. Default is 5.
}
  \item{numTops}{
	The number of pathway features used for feature selection. Default is 50.
}
  \item{iter}{
	The number of runs to split \code{xG} for optimal classifier selection. Default is 1.
}
  \item{Gamma}{
	A numeric value. The restart probability in directed random walk. Default is 0.7.
}
  \item{Alpha}{
	A proportional coefficient to balance the initial weights of genes and metabolites, which are used to construct the initial weights W0 for directed random walk.
}
  \item{fdr.output}{
	(Approximate) False Discovery Rate cutoff for output in significant genes table. Default is 0.2.
}
}
\details{
DRW-GM uses directed random walk to evaluate the topological importance of each gene in reconstructed gene-metabolite graph through integrating information from matched gene expression profiles and metabolomic profiles. The topological importance of genes are used to weight the genes for inferring robust DRW-GM-based pathway activities.
Then the pathway activities are selected to train the classifier.
}
\value{
	Fitted \code{"DRWPClassGM"} model object.
	\item{model}{Fitted \code{"fitModel"} model object.}
	\item{AUC}{The performance (AUC) of the classifier on \bold{feature selection set}.}
	\item{Accuracy}{The performance (Accuracy) of the classifier on \bold{feature selection set}.}
	\item{pathFeatures}{The selected pathway features to build the classifier.}
	\item{geneFeatures}{The genes used to infer the pathways in \code{pathFeatures}}
	\item{tScore}{The t statistic and p-value of each gene in \code{xG}}
	\item{vertexWeight}{The topological weights of vertexes in \code{globalGraph}}
	\item{pathSet}{The pathways used to construct the global directed gene-metabolite graph.}
	\item{globalGraph}{An \code{igraph} R object. The global directed gene-metabolite pathway graph.}
	\item{testStatistic}{The test method used to identify differential genes.}
	\item{classifier}{The method to train classifiers.}
	\item{nFolds}{The number of folds to split \code{xG}.}
	\item{numTops}{The number of pathway features used for feature selection.}
	\item{iter}{The number of runs to split \code{xG} for optimal classifier selection.}
	\item{Gamma}{The restart probability in directed random walk.}
	\item{Alpha}{The proportional coefficient to balance the initial weights of genes and metabolites.}
}
\references{
Liu, W., et al., Topologically inferring risk-active pathways toward precise cancer classification by directed random walk. 
  Bioinformatics, 2013. 29(17): p. 2169-77.
}
\author{
Wei Liu
}

\seealso{
\code{\link{predict.DRWPClassGM}}
}
\examples{
	data(GProf8511)
	data(MProf)
	data(pathSet)
	data(dGMGraph)
	fit <- fit.DRWPClassGM(xG=GProf8511$mRNA_matrix, yG.class1=GProf8511$normal, yG.class2=GProf8511$PCA,
                           xM=MProf$Meta_matrix, yM.class1=MProf$normal, yM.class2=MProf$PCA, DEBUG=TRUE, 
                           pathSet=pathSet, globalGraph=dGMGraph, testStatistic="t-test", classifier = "Logistic", 
                           normalize = TRUE, nFolds = 5, numTops=50, iter = 1, Gamma=0.7, Alpha = 0.5)
}