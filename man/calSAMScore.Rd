\name{calSAMScore}
\alias{calSAMScore}
\title{Calculate SAM score}
\description{Calculate the d statistic and q-value of each gene using the \code{samr} package.}
\usage{
calSAMScore(mRNA_matrix, normSample, diseaseSample, fdr.output = 0.05)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{fdr.output}{
    (Approximate) False Discovery Rate cutoff for output in significant genes table.
}
}
\details{
	The function calculates the d statistic and q-value of each gene in \code{mRNA_matrix} using the \code{samr} package. The output of the function will be used for topological weights and pathway activity inference.
}
\value{
	a p x 2 matrix is returned, which contains the d statistic and q-value of p genes in \code{mRNA_matrix}.
}
\references{
	Liu et al. Topologically inferring pathway activity toward precise cancer classification via integrating genomic and metabolomic data: prostate cancer as a case.
}
\author{
	Wei Liu <freelw@gmail.com>
}
\seealso{
	\code{\link{caltScore}}
}
\examples{
	# calculates the d statistic and q-value of each gene using SAM between normal samples and
	# PCA samples in GSE8511
	data(GProf8511)
	SAMScore <- calSAMScore(GProf8511$mRNA_matrix, GProf8511$normal, GProf8511$PCA, fdr.output = 0.2)
}
