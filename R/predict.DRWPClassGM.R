predict.DRWPClassGM <-
function(object, newx, type = c("class", "probability")){
	# browser()
	geneFeatures <- unique(unlist(object$geneFeatures))
	if(any(!(geneFeatures %in% rownames(newx)))){
		stop("The genes in newx do not match with that in the training set!")
	}
	newx.mu <- apply(newx,1,mean)
	newx.sigma <- apply(newx,1,sd)
	newx <- sweep(newx,1,newx.mu)
	newx <- sweep(newx,1,1/newx.sigma,'*')
	newx <- newx[geneFeatures,]
	pA <- getPathActivity(newx, object$pathSet[object$pathFeatures], object$vertexWeight, object$tScore)
	arff_newx <- data.frame(t(pA$pathActivity), check.names=F)
	predict(fit$model$model, newdata = arff_newx, type = type)
}
