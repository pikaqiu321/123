calSAMScore <-
function(mRNA_matrix, normSample, diseaseSample, fdr.output = 0.05){
	# calculate the t-test statistic of each gene in the expression profile
	# input:
	# mRNA_matrix: the expression profile
	# normSample: the index of normal samples
	# diseaseSample: the index of disease samples
	# output:
	# the t-test statistic of each gene in the expression profile
	# browser()
	tscore <- matrix(0, nrow = nrow(mRNA_matrix), ncol = 2)
	rownames(tscore) <- rownames(mRNA_matrix)
	colnames(tscore) <- c("tScore","p-value")
	tscore[,"p-value"] <- 1
	
	y <- rep(1,length(normSample)+length(diseaseSample))
	y[(length(normSample)+1):(length(normSample)+length(diseaseSample))] <- 2
	samfit<-SAM(mRNA_matrix,y,resp.type="Two class unpaired", fdr.output=fdr.output)
	
	for (i in 1 : samfit$siggenes.table$ngenes.up){
		idx <- which(rownames(mRNA_matrix) == samfit$siggenes.table$genes.up[i,"Gene Name"])
		tscore[idx,1] <- as.numeric(samfit$siggenes.table$genes.up[i,"Score(d)"])
		tscore[idx,2] <- as.numeric(samfit$siggenes.table$genes.up[i,"q-value(%)"])/100
	}
	
	for (i in 1 : samfit$siggenes.table$ngenes.lo){
		idx <- which(rownames(mRNA_matrix) == samfit$siggenes.table$genes.lo[i,"Gene Name"])
		tscore[idx,1] <- as.numeric(samfit$siggenes.table$genes.lo[i,"Score(d)"])
		tscore[idx,2] <- as.numeric(samfit$siggenes.table$genes.lo[i,"q-value(%)"])/100
	}	
	
	return(tscore)
}
