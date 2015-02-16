fit.DRWPClassGM <-
function(xG, yG.class1, yG.class2, xM=NULL, yM.class1=NULL, yM.class2=NULL, DEBUG=FALSE, pathSet, globalGraph, testStatistic=c("t-test","SAM"), classifier = "Logistic", normalize = TRUE, nFolds = 5, numTops=50, iter = 1, Gamma=0.7, Alpha = 0.5, fdr.output = 0.2){
	if(normalize){
		# normalize xG
		xG.mu <- apply(xG,1,mean)
		xG.sigma <- apply(xG,1,sd)	
		xG <- sweep(xG,1,xG.mu)
		xG <- sweep(xG,1,1/xG.sigma,'*')
		xG <- xG[,c(yG.class1, yG.class2)]
		yG.class1 <- 1:length(yG.class1)
		yG.class2 <-(length(yG.class1)+1):(length(yG.class1)+length(yG.class2))
		
		# normalize xM
		if(!is.null(xM)){
			xM.mu <- apply(xM,1,mean)
			xM.sigma <- apply(xM,1,sd)	
			xM <- sweep(xM,1,xM.mu)
			xM <- sweep(xM,1,1/xM.sigma,'*')
		}		
	}
	testStatistic <- match.arg(testStatistic)
	if(testStatistic == "t-test"){
		if(DEBUG) cat('Calculating t-test score...')	
		tScore <- caltScore(xG, yG.class1, yG.class2)
		if(DEBUG) cat('Done\n')
	}
	if(testStatistic == "SAM"){
		if(DEBUG) cat('Performing SAM...')	
		tScore <- calSAMScore(xG, yG.class1, yG.class2, fdr.output = fdr.output)	
		if(DEBUG) cat('Done\n')		
	}
	
	geneWeight <- -log(tScore[ ,2]+2.2e-16)
	geneWeight[which(is.na(geneWeight))] <- 0
	geneWeight <- (geneWeight - min(geneWeight)) / (max(geneWeight) - min(geneWeight))
	
	if(!is.null(xM)){
		metaWeight <- getMetaWeight(xM, yM.class1, yM.class2)
		metaWeight <- metaWeight[!is.na(metaWeight)]
		metaWeight <- (metaWeight - min(metaWeight)) / (max(metaWeight) - min(metaWeight))
		W0 <- getW0_GM(geneWeight=geneWeight * (1-Alpha), metaWeight=metaWeight * Alpha, globalGraph=globalGraph)
	}else{
		W0 <- getW0_GM(geneWeight=geneWeight, globalGraph=globalGraph)
	}
	
	if(DEBUG) cat('Performing directed random walk...')	
	vertexWeight <- DRW(igraphM = globalGraph, p0 = W0, EdgeWeight=FALSE, gamma = Gamma)
	names(vertexWeight) <- names(W0)
	if(DEBUG) cat('Done\n')	
	
	# gc()
	# library("RWeka")
	
	pA <- getPathActivity(xG, pathSet, vertexWeight, tScore)
	
	tScorePA <- caltScore(pA$pathActivity, yG.class1, yG.class2)
	tScore_pathway <- tScorePA[ ,1]
	Idx <- sort(tScore_pathway,decreasing = TRUE,index.return=TRUE)$ix
	tScore_pathway <- tScore_pathway[Idx]
	
	AUC <- 0
	Accuracy <- 0
	for(i in 1:iter){
		if(DEBUG) cat('iter..................................', i, '\n')
		for(j in 1:nFolds){
			if(DEBUG) cat('iter=', i,'Folds=', j, '\n')	
			yG.class1.test <- sample(yG.class1, floor(length(yG.class1) / nFolds))
			yG.class1.training <- setdiff(yG.class1, yG.class1.test)
			yG.class2.test <- sample(yG.class2, floor(length(yG.class2) / nFolds))
			yG.class2.training <- setdiff(yG.class2, yG.class2.test)
			
			pA.training <- pA$pathActivity[,c(yG.class1.training, yG.class2.training)]
			pA.test <- pA$pathActivity[,c(yG.class1.test, yG.class2.test)]
			
			yG.class1.training  <- 1:length(yG.class1.training)
			yG.class2.training <- (1+length(yG.class1.training)):(length(yG.class1.training)+length(yG.class2.training))
			yG.class1.test  <- 1:length(yG.class1.test)
			yG.class2.test <- (1+length(yG.class1.test)):(length(yG.class1.test)+length(yG.class2.test))
			
			classType_training <- rep(NA, dim(pA.training)[2])
			classType_training[yG.class1.training] <- "class1"
			classType_training[yG.class2.training] <- "class2"
			arffRW_training <- data.frame(t(pA.training), "class" = classType_training, check.names=F)
			
			classType_test <- rep(NA, dim(pA.test)[2])
			classType_test[yG.class1.test] <- "class1"
			classType_test[yG.class2.test] <- "class2"
			arffRW_test <- data.frame(t(pA.test), "class" = classType_test, check.names=F)
			
			# browser()
			resPredict <- fitModelGreedy(arffRW_training, arffRW_test, tScore_pathway, classifier = classifier, numTops=numTops)
			if(DEBUG){
				# print(i)
				print(list(model=resPredict$model, AUC=resPredict$AUC, pathFeatures = resPredict$pathFeatures, geneFeatures = pA$sigGenes[resPredict$pathFeatures]))
			}
			if(resPredict$AUC > AUC){
				# print(i)
				AUC <- resPredict$AUC
				Accuracy <- resPredict$Accuracy
				fit <- list(model=resPredict$model, AUC=resPredict$AUC, Accuracy=resPredict$Accuracy, pathFeatures = resPredict$pathFeatures, geneFeatures = pA$sigGenes[resPredict$pathFeatures], tScore = tScore, vertexWeight = vertexWeight, pathSet = pathSet, globalGraph = globalGraph, testStatistic = testStatistic, classifier = classifier, nFolds = nFolds, numTops=numTops, iter = iter, Gamma=Gamma, Alpha = Alpha)
			}else{
				if(resPredict$AUC == AUC & resPredict$Accuracy > Accuracy){
					# print(i)
					AUC <- resPredict$AUC
					Accuracy <- resPredict$Accuracy
					fit <- list(model=resPredict$model, AUC=resPredict$AUC, Accuracy=resPredict$Accuracy, pathFeatures = resPredict$pathFeatures, geneFeatures = pA$sigGenes[resPredict$pathFeatures], tScore = tScore, vertexWeight = vertexWeight, pathSet = pathSet, globalGraph = globalGraph, testStatistic = testStatistic, classifier = classifier, nFolds = nFolds, numTops=numTops, iter = iter, Gamma=Gamma, Alpha = Alpha)
				}
			}
		}
	}
	class(fit) <- "DRWPClassGM"
	return(fit)
}
