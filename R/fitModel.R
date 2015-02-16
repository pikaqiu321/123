fitModel <-
function(expr_training, expr_test, classifier){
	# library("RWeka")
	switch(classifier, "libsvm" = {
		WPM("load-package", "LibSVM")
		libsvm <- make_Weka_classifier("weka/classifiers/functions/LibSVM")
		# print("libsvm Supplied test set")
		model <- libsvm(class~., data = expr_training)
		eTestSet <- evaluate_Weka_classifier(model, newdata = expr_test, class = TRUE)
		# print(eTestSet)
	},"libLINEAR" = {
		WPM("load-package", "LibLINEAR")
		libLINEAR <- make_Weka_classifier("weka/classifiers/functions/LibLINEAR")
		# print("libLINEAR Supplied test set")
		model <- libLINEAR(class~., data = expr_training)
		eTestSet <- evaluate_Weka_classifier(model, newdata = expr_test, class = TRUE)
		# print(eTestSet)
	},"NaiveBayes" = {
		# WPM("load-package", "LibSVM")
		# libsvm <- make_Weka_classifier("weka/classifiers/functions/LibSVM")
		NaiveBayes <- make_Weka_classifier("weka/classifiers/bayes/NaiveBayes")
		# print("NaiveBayes Supplied test set")
		model <- NaiveBayes(class~., data = expr_training)
		eTestSet <- evaluate_Weka_classifier(model, newdata = expr_test, class = TRUE)
		# # print(eTestSet)
	}, 'J48' = 	{
		# print("J48 Supplied test set")
		model <- J48(class~., data = expr_training)
		eTestSet <- evaluate_Weka_classifier(model, newdata = expr_test, class = TRUE)
		# # print(eTestSet)
		# return (eTestSet)		
	}, 'Logistic' = {
		# print("Logistic Supplied test set")
		model <- Logistic(class~., data = expr_training)
		eTestSet <- evaluate_Weka_classifier(model, newdata = expr_test, class = TRUE)	
		# print(eTestSet)
		# browser()
		# return (eTestSet)
	})
	fit <- list(model=model, eTestSet=eTestSet)
	class(fit) <- "fitModel"
	return (fit)	
}
