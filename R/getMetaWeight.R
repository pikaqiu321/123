getMetaWeight <-
function(metaProf, normMetaSmpl, diseaseMetaSmpl){
	# get the initial weight of metabolites
	# metaProf: metabolite profile
	# normMetaSmpl: index of normal samples
	# diseaseMetaSmpl: index of disease samples
	metaProf <- as.matrix(metaProf)
	pValue <- rep(NA, dim(metaProf)[1])
	names(pValue) <- rownames(metaProf)
	for (i in 1 : dim(metaProf)[1]){
		tt <- wilcox.test(metaProf[i,normMetaSmpl],metaProf[i,diseaseMetaSmpl])
		pValue[i] <- tt$p.value
	}
	
	metaWeight <- -log(pValue)
	
	metaWeight[which(pValue < 0.01)] <- 1
	metaWeight[which(pValue >= 0.01)] <- 0
	
	return(metaWeight)	
}
