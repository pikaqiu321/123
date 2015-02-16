evaluate.DRWPClassGM <-
function(object, newx, newy.class1, newy.class2){
	geneFeatures <- unique(unlist(object$geneFeatures))
	if(any(!(geneFeatures %in% rownames(newx)))){
		stop("The genes in newx do not match with those in the training set!")
	}
	newx.mu <- apply(newx,1,mean)
	newx.sigma <- apply(newx,1,sd)
	newx <- sweep(newx,1,newx.mu)
	newx <- sweep(newx,1,1/newx.sigma,'*')
	newx <- newx[geneFeatures,c(newy.class1,newy.class2)]
	newy.class1 <- 1:length(newy.class1)
	newy.class2 <- (length(newy.class1)+1):(length(newy.class1)+length(newy.class2))
	
	pA <- getPathActivity(newx, object$pathSet[object$pathFeatures], object$vertexWeight, object$tScore)
	
	classType_test <- rep(NA, dim(newx)[2])
	classType_test[newy.class1] <- "class1"
	classType_test[newy.class2] <- "class2"
	arff_newx <- data.frame(t(pA$pathActivity), "class" = classType_test, check.names=F)
	evaluate_Weka_classifier(object$model$model, newdata = arff_newx, class = TRUE)
}
