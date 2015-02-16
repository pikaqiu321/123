fitModelGreedy <-
function(expr_training, expr_test, TValue, classifier = "Logistic", numTops = 50){
	res.fitModel <- fitModel(expr_training[, c(names(TValue)[1],"class")], expr_test[, c(names(TValue)[1],"class")], classifier)
	eString <- unlist(strsplit(res.fitModel$eTestSet$string,"\n"))
	eString.4 <- unlist(strsplit(eString[4], "[[:blank:]]+"))
	Accuracy <- as.numeric(eString.4[5])
	eString.20 <- unlist(strsplit(eString[20],"[[:blank:]]+"))
	AUC <- as.numeric(eString.20[9])
	flag <- 1
	iter <- min(numTops, ncol(expr_training) - 1)
	# print(ncol(expr_training) - 1)
	for (i in 2 : iter)	{
		res.fitModel.Tmp <- fitModel(expr_training[, c(names(TValue)[c(flag, i)],"class")], expr_test[, c(names(TValue)[c(flag, i)],"class")], classifier)
		eStringTmp <- unlist(strsplit(res.fitModel.Tmp$eTestSet$string,"\n"))
		eStringTmp.4 <- unlist(strsplit(eStringTmp[4], "[[:blank:]]+"))
	    AccuracyTmp <- as.numeric(eStringTmp.4[5])
		eStringTmp.20 <- unlist(strsplit(eStringTmp[20],"[[:blank:]]+"))
		AUCTmp <- as.numeric(eStringTmp.20[9])
		if (AUC < AUCTmp){
			flag <- c(flag, i)
			res.fitModel <- res.fitModel.Tmp
			Accuracy <- AccuracyTmp
			AUC <- AUCTmp
		}
		if(AUC == AUCTmp & Accuracy < AccuracyTmp){
			flag <- c(flag, i)
			res.fitModel <- res.fitModel.Tmp
			Accuracy <- AccuracyTmp
			AUC <- AUCTmp
		}
	}
	return(list(model=res.fitModel, AUC=AUC, Accuracy=Accuracy, pathFeatures=names(TValue)[flag]))
}
