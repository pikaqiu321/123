Package DRWPClassGM implements the DRW-GM method. DRW-GM is a pathway-based classification method which performs classifier construction and precise disease status prediction by joint analysis of genomic and metabolomic data and pathway topology.
How to use DRWPClassGM
Essential dependence platform
To run DRWPClassGM, you should make sure that your computer has installed Weka. The Weka downloaded website (in our analysis we installed weka-3-7-12jre.exe):
http://www.cs.waikato.ac.nz/~ml/weka/downloading.html

The available resource of DRWPClass can be downloaded as follows:
Package source:	DRWPClassGM_1.0.tar.gz 

Windows binary:	DRWPClassGM_1.0.zip 

Reference manual:	DRWPClassGM-manual.pdf 

http://222.170.78.208/DRW-GM/

Dependent packages:
Depends: R (>= 3.1.1), igraph, Matrix, RWeka, samr
###########Example code for predicting disease class##############
library(DRWPClassGM)
# load example data
data(GProf8511)
data(GProf3325)
data(MProf)
data(pathSet)
data(dGMGraph)
# train a classifier
fit <- fit.DRWPClassGM(xG=GProf8511$mRNA_matrix, yG.class1=GProf8511$normal, yG.class2=GProf8511$PCA, xM=MProf$Meta_matrix, yM.class1=MProf$normal, yM.class2=MProf$PCA, DEBUG=TRUE, pathSet=pathSet, globalGraph=dGMGraph, testStatistic="t-test", classifier = "Logistic")
# prediction
predict.DRWPClassGM(object=fit, newx=GProf3325$mRNA_matrix[,c(GProf3325$normal,GProf3325$PCA)], type = "class")
# evaluate classification performance
evaluate.DRWPClassGM(object=fit, newx=GProf3325$mRNA_matrix, newy.class1=GProf3325$normal, newy.class2=GProf3325$PCA)
