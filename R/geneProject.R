geneProject = function(train.data, test.data, gene.sign=NULL, normalize=c("zscore", "combat", "none"), PC=1, return.data=FALSE) {
  
  ## Takes in gene expression datasets with genes as rows and samples as columns
  ## and builds a model to score new samples using SVD
  ## Written by: Josh Campbell, Kahkeshan Hijazi
  ## 1-15-2012
  
  ## Parameters: #######################################################################
  #
  # train.data		Matrix or data frame. Genes as rows and samples as columns.
  # test.data		Matrix or data frame. New samples to score.
  # gene.sign		Numeric vector or NULL. The sign of the gene with respect to the phenotype.
  #					Positive values indicate that the gene was positively correlated with
  #					the phenotype (i.e. up-regulated).  Negative values indicate the gene
  #					was negatively correlated with the phenotype (i.e. down-regulated).  These
  #					values will be used to try and give a direction to the PCs.  If NULL, the
  #					sign of the rotation matrix will not be altered.  Values could be  
  #					t-statistics or correlation coefficients for each gene. Could also be a
  #					1 for "up" genes and -1 for "down" genes.  Genes with values that are
  #					NA, NULL, Inf, or 0 will be filtered out after normalization.
  # normalize		One of zscore, combat, or none.  Note that if combat is used, all
  #					overlapping genes in the two datasets should be used as input.
  # PC				Integer or vector of integers.  Which PCs to return as the sample score.
  # return.data  	Boolean.  Whether or not to return the normalized data used in the 
  #					building of the model.
  #
  ######################################################################################
  
  ## Values: ###########################################################################
  #
  # svd					The object returned from the "prcomp" function 
  # signs.reversed		Boolean indicating whether the signs of the "rotation" matrix were
  #						reversed based on the gene.sign vector
  # train.scores		Predicted scores of the train samples
  # test.scores			Predicted scores of the test samples
  # svd.gene.sign		Vector or 1 or -1 indicating if each gene was apart of the positive
  #						or negative group used to orient the PC socres
  # normalize.method	Normalization method used.
  # train.normalized	Normalized training data actually used to build the model
  # test.normalized		Normalized test data where the scores were predicted  
  #
  ######################################################################################
  
  normalize = match.arg(normalize)
  
  ## Check to make sure test and train datasets are the same size
  train.nrow = nrow(train.data)
  test.nrow = nrow(test.data)
  if(train.nrow != test.nrow) {
    stop("Test dataset must have the same number of genes as the training dataset")
  }
  
  ## Check to see if row IDs are the same across datasets.  prcomp will throw an
  ## error if the IDs do not match up when predicting in a new group
  if(sum(rownames(train.data) == rownames(test.data)) != train.nrow) {
    warning("The test dataset does not have the same row ids as the training dataset. Setting row ids for the test data to those of the training data")
	rownames(test.data) = rownames(test.data)
  }

  ## Process gene.sign variable
  if(is.null(gene.sign)) {
    all.gene = 1:nrow(train.data)
	new.gene.sign = NULL
  } else {
    gene.sign = as.numeric(gene.sign)
    pos.ind = which(gene.sign > 0)
    neg.ind = which(gene.sign < 0)
    null.ind = which(is.null(gene.sign) | is.na(gene.sign) | is.infinite(gene.sign) | gene.sign == 0)
    all.gene = sort(c(pos.ind, neg.ind))
	new.gene.sign = sign(gene.sign[all.gene])
  }
  
  ## Normalize data so projection can be done across datasets
  if(normalize == "zscore") {
    
	train.data = scale(t(train.data[all.gene,]))
	test.data = scale(t(test.data[all.gene,]))
	svd.center=FALSE
	svd.scale=FALSE
	
  } else if (normalize == "combat") {
    
	require(sva)
	batch = as.factor(rep(c("train","test"), c(ncol(train.data), ncol(test.data))))
	mod0 = rep(1, length(batch))
	all.data = cbind(train.data, test.data)
	cat("Performing batch correction using ComBat...\n")
	norm.data = ComBat(dat=all.data, batch=batch, mod=mod0, numCovs=NULL, par.prior=TRUE)
	train.data = t(norm.data[all.gene,batch=="train"])
	test.data = t(norm.data[all.gene,batch=="test"])
	rm(all.data)
	svd.center=TRUE
	svd.scale=TRUE
	
  } else {
    
	train.data = t(train.data[all.gene,])
	test.data = t(test.data[all.gene,])
	svd.center=FALSE
	svd.scale=FALSE
}
  
  ## Perform SVD.  
  svd = prcomp(train.data, center=svd.center, scale=svd.scale)
  
  ## Set these to FALSE just in case.  Found a bug where sometimes
  ## the function would fill in these parameters even though they
  ## were set to FALSE in the function call (which would mess up
  ## the prediction)
  if(svd.center == FALSE) {
    svd$center=FALSE   
  }
  if(svd.scale == FALSE) {
    svd$scale=FALSE
  }
  
  ## Examines the signs of the first PC to see if the positive genes
  ## are indeed positive in the model.  Same for negative genes.
  ## If not, then the signs of the rotation matrix are reversed.
  signs.reversed = FALSE
  if(!is.null(new.gene.sign)) {
    pos.mean = mean(svd$rotation[new.gene.sign == 1,1])
    neg.mean = mean(svd$rotation[new.gene.sign == -1,1])
	
	if(pos.mean < 0 & neg.mean > 0) {
	  svd$rotation = -svd$rotation
	  signs.reversed=TRUE
	}
  }
  
  ## Make predictions
  pred.train = predict(svd, train.data)
  pred.test = predict(svd, test.data)
  
  if(return.data == TRUE) {
    return(list("svd"=svd, "train.scores"=pred.train[,PC], "test.scores"=pred.test[,PC], "signs.reversed"=signs.reversed, "normalize.method"=normalize, "svd.gene.sign"=new.gene.sign, "train.normalized"=t(train.data), "test.normalized"=t(test.data)))
  } else {
    return(list("svd"=svd, "train.scores"=pred.train[,PC], "test.scores"=pred.test[,PC], "signs.reversed"=signs.reversed, "normalize.method"=normalize, "svd.gene.sign"=new.gene.sign))
  }
}
