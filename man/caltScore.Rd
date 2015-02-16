\name{caltScore}
\alias{caltScore}
\title{	Calculate tScore.}
\description{
	Calculate the t statistic and p-value of each gene in \code{mRNA_matrix} using the \code{t-test}.
}
\usage{
caltScore(mRNA_matrix, normSample, diseaseSample)
}
\arguments{
  \item{mRNA_matrix}{
	a p x n matrix of expression measurements with p genes and n samples.
}
  \item{normSample}{
	a integer vector comprising the indexes of the normal (first class) samples in \code{mRNA_matrix}.
}
  \item{diseaseSample}{
	a integer vector comprising the indexes of the disease (second class) samples in \code{mRNA_matrix}.
}
}
\details{
	The function calculates the t statistic and p-value of each gene in \code{mRNA_matrix} using \code{t-test}. The output of the function will be used for topological weights and pathway activity inference.
}
\value{
	a p x 2 matrix is returned, which contains the t statistic and p-value of p genes in \code{mRNA_matrix}.
}
\references{
	Liu et al. Topologically inferring pathway activity toward precise cancer classification via integrating genomic and metabolomic data: prostate cancer as a case.
}
\author{
	Wei Liu <freelw@gmail.com>
}
\seealso{
	\code{\link{calSAMScore}}
}
\examples{
	# calculates the t statistic and p-value of each gene using t-test between normal samples and
	# PCA samples in GSE8511
	data(GProf8511)
	tScore <- caltScore(GProf8511$mRNA_matrix, GProf8511$normal, GProf8511$PCA)
}