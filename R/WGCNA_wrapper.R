
plotSoftThreshold = function(datExpr, powerVector = 1:20, RsquaredCut = 0.90, corFnc = "cor", networkType = "unsigned", savePlots=NULL) {
  
  ## Tries to automatically identify the best softPower to use in the WGCNA pipeline
  ## This function will make several plots which can be helpful in choosing the softPower
  ## and will try to make a guess at the best one.  The plots should be examined to see
  ## if the guess looks good.
  ##
  ## Written by Josh Campbell
  ## 1-15-2013
  
  ## Parameters: #######################################################################
  #
  # datExpr 			Matrix or data frame. Genes as rows and samples as columns.
  # 
  # For powerVector, corFnc, networkType, and RsquaredCut, see functions adjacency, 
  # pickSoftThreshold, softConnectivity, and scaleFreePlot in the WGCNA package and
  # tutorials for more info.
  #
  ######################################################################################
  
  ## Values: ###########################################################################
  #
  # powerEstimate		Best guess at the softPower parameter
  #
  ######################################################################################
  
  require(WGCNA)  
  datExpr = t(datExpr)
  options(stringsAsFactors=FALSE)
  
  # Call the network topology analysis function
  sft = pickSoftThreshold(datExpr, powerVector = powerVector, RsquaredCut = RsquaredCut, corFnc = corFnc, networkType = networkType, verbose = 5)

  if(!is.null(savePlots)) {
    pdf(savePlots, useDingbats=FALSE)
    cex1 = 0.9;
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powerVector,cex=cex1,col="red");

    # this line corresponds to using an R^2 cut-off of h
    abline(h=0.90,col="red")

    # Mean connectivity as a function of the soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powerVector, cex=cex1,col="red")
  
    for(i in powerVector) {
      softConn = softConnectivity(datExpr, corFnc = corFnc, type = networkType, power=i)
	  scaleFreePlot(softConn, main=paste("Power=", i, sep=""))
    }
    dev.off()
  }
  return(sft$powerEstimate)
}





wgcna = function(datExpr, softPower = 6, corFnc="cor", networkType="unsigned", TOMType="unsigned", cutHeight=0.995, minModuleSize=20, MEDissThres=0.15, saveAdjacency=NULL, saveDissTOM=NULL, plotModuleHeatmaps=NULL) {
  
  ## Wrapper for the WGCNA function which can be used to identify coeexpression modules within a dataset
  ## Many different parameters can be tweeked for WGCNA.  These parameters are the ones used by
  ## Bin Zhang who helped write WGCNA and is now at Eric Schadt's group at Mt. Sinai.  They are the same
  ## default parameters used in the blockwiseModules
  ##
  ## For more information on WGCNA, how it works and what it's paramters mean, see:
  ## http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/
  ## and 
  ## http://www.biomedcentral.com/1471-2105/9/559
  ##  
  ## Written by Josh Campbell
  ## 1-15-2013
  
  ## Parameters: #######################################################################
  #
  # datExpr 			Matrix or data frame. Genes as rows and samples as columns.
  # saveAdjacency		NULL or string.  If string, the adajacency matrix will be saved
  #						into an R object with name saveAdjacency
  # saveDissTOM			NULL or string.  If string, the TOM dissimilatrity matrix will
  #						be saved into an R object with name saveDissTOM
  # plotModuleHeatmaps	NULL or string.  If string, a pdf with the name plotModuleHeatmaps
  #						will be made with a heatmap for each module.  Requires heatmap3
  #						and gplots.  
  #
  # For softPower, corFnc, networkType, TOMType, cutHeight, minModuleSize,
  # and MEDissThres, see functions adjacency, TOMdist, cutreeDynamic, and
  # mergeCloseModules in the WGCNA package and tutorials for more info.
  #
  ######################################################################################
  
  ## Values: ###########################################################################
  #
  # modules				A data frame with one row for each gene giving the gene ID and
  #						which module each gene belongs to.  Note that module "grey"
  #						with ID 0 is a "junk" module for all genes that could not
  #						be clustered.  
  # MEs					Matrix with one row for each sample and one columne for each module.
  #						Values are the "eigen gene" for that module in that sample.
  #
  ######################################################################################
  
  require(WGCNA)
  options(stringsAsFactors=FALSE)
  
  cat("Calculating adjacency matrix...\n")
  datExpr = t(datExpr)
  adjacency = adjacency(datExpr, power = softPower, corFnc=corFnc, type=networkType);

  
  cat("Calculating TOM matrix...\n")
  dissTOM = TOMdist(adjacency, TOMType=TOMType);
  
  if(is.character(saveAdjacency)) {
    save(adjacency, file=saveAdjacency)
  }
  rm(adjacency)
  
  cat("Performing garbage collection...\n")
  gc()

  if(is.character(saveDissTOM)) {
    save(dissTOM, file=saveDissTOM)
  }  

  cat("Performing garbage collection...\n")
  gc()

  cat("Clustering genes using TOM matrix...\n")
  geneTree = flashClust(as.dist(dissTOM), method = "average");

  dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight=0.995, method="tree", minClusterSize = minModuleSize);
  dynamicColors = labels2colors(dynamicMods)

  cat("Merging close modules...\n")
  merged.modules = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 5)
  mergedColors = merged.modules$colors;
  mergedMEs = merged.modules$newMEs;

  colorOrder = names(sort(table(mergedColors), decreasing=TRUE));
  colorOrder = c("grey", setdiff(colorOrder, "grey"))
  moduleIDs = match(mergedColors, colorOrder)-1;
  
  if(is.character(plotModuleHeatmaps)) {
    require(heatmap3)
	require(gplots)
	cat("Creating module heatmaps...\n")
	
	pdf(plotModuleHeatmaps)
	for(i in setdiff(unique(sort(moduleIDs)), 0)) {
	  
	  index = which(moduleIDs == i)
	  current.color = unique(mergedColors[index])
	  title = paste("Module ID: ", i, "; Module Color: ", current.color, sep="")
	  
	  me = paste("ME", current.color, sep="")
	  heatmap3(t(datExpr[order(mergedMEs[,me]),index]), col=bluered(100), row.dendrogram=F, Colv=NA, labRow=F, labCol=F, main=title)
	}
	dev.off()
  }
  
  cat("Completed!\n")
  return(list(modules=data.frame("Gene"=colnames(datExpr), "ID"=moduleIDs, "Color"=mergedColors), "MEs"=mergedMEs))
}












wgcnaReturnTOM = function(datExpr, softPower = 6, corFnc="cor", networkType="unsigned", TOMType="unsigned", cutHeight=0.995, minModuleSize=20, MEDissThres=0.15) {
  
  ## Wrapper for the WGCNA function which can be used to identify coeexpression modules within a dataset
  ## Many different parameters can be tweeked for WGCNA.  These parameters are the ones used by
  ## Bin Zhang who helped write WGCNA and is now at Eric Schadt's group at Mt. Sinai.  They are the same
  ## default parameters used in the blockwiseModules
  ##
  ## For more information on WGCNA, how it works and what it's paramters mean, see:
  ## http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/
  ## and 
  ## http://www.biomedcentral.com/1471-2105/9/559
  ##  
  ## Written by Josh Campbell
  ## 1-15-2013
  
  ## Parameters: #######################################################################
  #
  # datExpr 			Matrix or data frame. Genes as rows and samples as columns.
  # saveAdjacency		NULL or string.  If string, the adajacency matrix will be saved
  #						into an R object with name saveAdjacency
  # saveDissTOM			NULL or string.  If string, the TOM dissimilatrity matrix will
  #						be saved into an R object with name saveDissTOM
  # plotModuleHeatmaps	NULL or string.  If string, a pdf with the name plotModuleHeatmaps
  #						will be made with a heatmap for each module.  Requires heatmap3
  #						and gplots.  
  #
  # For softPower, corFnc, networkType, TOMType, cutHeight, minModuleSize,
  # and MEDissThres, see functions adjacency, TOMdist, cutreeDynamic, and
  # mergeCloseModules in the WGCNA package and tutorials for more info.
  #
  ######################################################################################
  
  ## Values: ###########################################################################
  #
  # modules				A data frame with one row for each gene giving the gene ID and
  #						which module each gene belongs to.  Note that module "grey"
  #						with ID 0 is a "junk" module for all genes that could not
  #						be clustered.  
  # MEs					Matrix with one row for each sample and one columne for each module.
  #						Values are the "eigen gene" for that module in that sample.
  #
  ######################################################################################
  
  require(WGCNA)
  options(stringsAsFactors=FALSE)
  
  cat("Calculating adjacency matrix...\n")
  datExpr = t(datExpr)
  adjacency = adjacency(datExpr, power = softPower, corFnc=corFnc, type=networkType);

  
  cat("Calculating TOM matrix...\n")
  dissTOM = TOMdist(adjacency, TOMType=TOMType);
  
  rm(adjacency)
  
  cat("Performing garbage collection...\n")
  gc()

  cat("Performing garbage collection...\n")
  gc()

  cat("Clustering genes using TOM matrix...\n")
  geneTree = flashClust(as.dist(dissTOM), method = "average");

  dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight=0.995, method="tree", minClusterSize = minModuleSize);
  dynamicColors = labels2colors(dynamicMods)

  cat("Merging close modules...\n")
  merged.modules = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 5)
  mergedColors = merged.modules$colors;
  mergedMEs = merged.modules$newMEs;

  colorOrder = names(sort(table(mergedColors), decreasing=TRUE));
  colorOrder = c("grey", setdiff(colorOrder, "grey"))
  moduleIDs = match(mergedColors, colorOrder)-1;
  
  cat("Completed!\n")
  return(list(modules=data.frame("Gene"=colnames(datExpr), "ID"=moduleIDs, "Color"=mergedColors), "MEs"=mergedMEs, "TOM"=1-dissTOM))
}












wgcnaModuleDC = function(datExpr1, datExpr2, powerVector = 1:20, RsquaredCut = 0.90, corFnc = "cor", networkType = "unsigned", TOMType="unsigned", cutHeight=0.995, minModuleSize=20, MEDissThres=0.15) {
  datExpr1 = as.matrix(t(datExpr1))
  datExpr2 = as.matrix(t(datExpr2))
  
  ## Identify modules
  datExpr1.power = plotSoftThreshold(datExpr1, powerVector=powerVector, RsquaredCut=RsquaredCut, corFnc=corFnc, networkType=networkType)
  datExpr2.power = plotSoftThreshold(datExpr2, powerVector=powerVector, RsquaredCut=RsquaredCut, corFnc=corFnc, networkType=networkType)
  
  datExpr1.wgcnaWithTOM = wgcna(datExpr1, power=datExpr1.power, corFnc=corFnc, networkType=networkType, TOMType=TOMType, cutHeight=cutHeight, minModuleSize=minModuleSize, MEDissThres=MEDissThres)
  datExpr2.wgcnaWithTOM = wgcna(datExpr2, power=datExpr2.power, corFnc=corFnc, networkType=networkType, TOMType=TOMType, cutHeight=cutHeight, minModuleSize=minModuleSize, MEDissThres=MEDissThres)
  
  ## Calculate real moduleDC stat
  stat1 = moduleDCstatistic(dissTOM1, dissTOM2, modules1)
  stat2 = moduleDCstatistic(dissTOM2, dissTOM1, modules2)
  
  ## Run permuation using WGCNA
      


}



clusterMEs = function(MEs, coln) {
  require(mclust)
  clust.num = c()
  clust.labels = c()
  tert.labels = c()
  all.models = list()
  for(i in 1:ncol(MEs)) {
    model = Mclust(MEs[,i], verbose = FALSE)
    clust.num = c(clust.num, model$G)
    clust.labels = rbind(clust.labels, model$classification)
    all.models = c(all.models, list(model))
    
    tert = as.numeric(cut(MEs[,i], quantile(MEs[,i], c(0, 0.33, 0.66, 1)), include.lowest = TRUE))
    tert.labels = rbind(tert.labels, tert)
  }
  rownames(clust.labels) = colnames(MEs)
  colnames(clust.labels) = coln
  rownames(tert.labels) = colnames(MEs)
  colnames(tert.labels) = coln
  return(list(tert=tert.labels, clust=clust.labels))
}


fill = function(v, n) { return(c(v, rep("", n-length(v))))}

createWGCNATable = function(res, norm, entrez, symbol, prefix, n=nrow(norm)) {
  
  summary.symbol = c()
  summary.entrez = c()
  for(i in colnames(res$MEs)) {
    ix.color = gsub("ME", "", i)
    ix = res$modules$Color == ix.color  
    
    temp.cor = cor(res$MEs[,i], t(norm[ix,]))
    direction = sign(temp.cor)
    
    pos.ix = which(ix)[direction == 1]
    neg.ix = which(ix)[direction == -1]
    pos.ix = pos.ix[order(temp.cor[direction == 1], decreasing=T)]
    neg.ix = neg.ix[order(temp.cor[direction == -1], decreasing=F)]
    
    pos = setdiff(symbol[pos.ix], "")
    neg = setdiff(symbol[neg.ix], "")
    temp = cbind(fill(pos, nrow(norm)), fill(neg, nrow(norm)))
    colnames(temp) = paste(ix.color, c("positive", "negative"), sep="_")
    summary.symbol = cbind(summary.symbol, temp)
    
    pos = setdiff(entrez[pos.ix], "")
    neg = setdiff(entrez[neg.ix], "")
    temp = cbind(fill(pos, nrow(norm)), fill(neg, nrow(norm)))
    colnames(temp) = paste(ix.color, c("positive", "negative"), sep="_")
    summary.entrez = cbind(summary.entrez, temp)
    
  }  
  write.table(summary.symbol, paste0(prefix, "_ME_Genes_Symbol.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(summary.entrez, paste0(prefix, "_ME_Genes_Entrez.txt"), sep="\t", row.names=FALSE, quote=FALSE)
  
  return(list(symbol=summary.symbol, entrez=summary.entrez))
}


