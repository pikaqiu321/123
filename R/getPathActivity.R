getPathActivity <-
function(x, pathSet, w, vertexZP){
	# infer pathway expression profile
    # gNonMetabolic: nonMetabolic pathway
	# gMetabolic: metabolic pathway
	# mRNA_matrix_training: the expression profile of training set
	# mRNA_matrix_test: the expression profile of feature evaluation set
	# mRNA_matrix_crosstest: the expression profile of test set
	# vertexWeight: the weight of genes from random walk
	# vertexTScore: the t-test statistic of each gene in the global pathway network
	# output:
	# the pathway expression profile
	
	# infer non metabolic pathway expression profile
	pathActivity <- c()
	sigGenes <- vector("list", length = length(pathSet))  # The differential genes to construct pathActivity
	names(sigGenes) <- names(pathSet)
	TValue.pathActivity <- matrix(NA, nrow = length(pathSet), ncol = 1)
	rownames(TValue.pathActivity) <- names(pathSet)	
	for (i in 1 : length(pathSet)){
		Vpathwayi <- pathSet[[i]]
		if (length(Vpathwayi) > 0){
			n <- 0    # the number of differential genes in ith pathway 			
			pathActivity_tmp <- matrix(nrow=1,ncol=dim(x)[2],data=0) 
			TValue.pathActivity_tmp <- 0
			sigGenesi <- c()
			Idx_pathwayi <- c()   
			for (j in 1 : length(Vpathwayi)){
				Idx <- which(rownames(x)==Vpathwayi[j])
				if (length(Idx) > 0){
					if ( rownames(x)[Idx] %in% names(w)){
						if(vertexZP[rownames(x)[Idx],"p-value"] < 0.05){
							pathActivity_tmp <- pathActivity_tmp + sign(vertexZP[rownames(x)[Idx],"tScore"]) * w[rownames(x)[Idx]] * x[Idx,]							
							n <- n + 1
							Idx_pathwayi <- rbind(Idx_pathwayi,Idx)	
							sigGenesi <- c(sigGenesi, Vpathwayi[j])
						}
					}
				}
			}
			if(n > 0){
				pathActivity_tmp <- pathActivity_tmp / sqrt(sum(w[rownames(x)[Idx_pathwayi]]^2))
				rownames(pathActivity_tmp) <- names(pathSet)[i] 
				pathActivity <- rbind(pathActivity, pathActivity_tmp)
				sigGenes[[i]] <- sigGenesi
			}#else
		}#else
	}
	colnames(pathActivity) <- colnames(x)
	return(list(pathActivity=pathActivity, sigGenes=sigGenes))	
}
