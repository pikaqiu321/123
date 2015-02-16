getW0_GM <-
function(geneWeight, metaWeight=NULL, globalGraph){
	Vertexs <- V(globalGraph)
	W0 <- rep(0, length(Vertexs))
	names(W0) <- Vertexs$name
	if(!is.null(metaWeight)){
		for(i in 1:length(metaWeight)){
			idx <- grep(names(metaWeight[i]),names(W0))
			if(length(idx) > 0){
				W0[idx] <- metaWeight[i]
			}	
		}
	}
	for(i in 1 : length(W0)){
		idx <- which(names(geneWeight) == names(W0[i]))
		if(length(idx) > 0){
			W0[i] <- geneWeight[idx]
		}
	}
	# browser()
	W0 <- W0/sum(W0)
}
