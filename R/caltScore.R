caltScore <-
function(mRNA_matrix, normSample, diseaseSample){
	# calculate the t-test statistic of each gene in the expression profile
	# input:
	# mRNA_matrix: the expression profile
	# normSample: the index of normal samples
	# diseaseSample: the index of disease samples
	# output:
	# the t-test statistic of each gene in the expression profile
	# browser()
	tscore <- matrix(NA, nrow = nrow(mRNA_matrix), ncol = 2)
	rownames(tscore) <- rownames(mRNA_matrix)
	colnames(tscore) <- c("tScore","p-value")
	for ( i in 1 : nrow(mRNA_matrix)){
	    tscore_tmp <- t.test(mRNA_matrix[i,diseaseSample], mRNA_matrix[i,normSample], var.equal=TRUE)
		tscore[i,1] <- tscore_tmp$statistic
		tscore[i,2] <- tscore_tmp$p.value
	}
	return(tscore)
}
