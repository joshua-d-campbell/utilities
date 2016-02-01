maf2tabler = function(maf, variant.categories=tcga.mapping(), sample.list=NULL, gene.list=NULL, gene.order="frequency", sample.order="mutation.type", sample.column="Tumor_Sample_Barcode", variant.column="Variant_Classification", gene.column="Hugo_Symbol", na.samples=NULL, na.genes=NULL, verbose=TRUE) {
 
  ## Pull relevant columns from maf file
  samples = as.character(maf[,sample.column])
  genes = as.character(maf[,gene.column])
  variants = as.character(maf[,variant.column])
  
  ## Set up sample list
  s1 =  setdiff(sample.list, samples)
  s2 =  setdiff(samples, sample.list)
  if(!is.null(sample.list) & length(s1) > 0) {
    if(verbose == TRUE) { cat("The following samples were not in the maf:", "\n", paste(s1, collapse="\n"), "\n\nThey will be given zero counts across all genes\n") }
    samples.list = unique(c(sample.list, samples))
  } else if (!is.null(sample.list) & length(s2) > 0) {
    if(verbose == TRUE) { cat("The following samples in the maf were not in sample.list:", "\n", paste(s2, collapse="\n"), "\n\nThey will be excluded from the matrix and variant counts.\n") }
     maf = subset(maf, maf[,sample.column] %in% sample.list)
  } else {
    sample.list = unique(samples)
  }
  
  ## Set up gene list
  if(!is.null(gene.list)) {
    g1 = setdiff(gene.list, genes)
    if(length(g1) > 0) {
      if(verbose == TRUE) { cat("The following genes were not in the maf:", "\n", paste(genes.extra, collapse="\n"), "\n") }
    }
  } else {
    gene.list = unique(genes)
  }

 
  ## Generate summary counts for mutations by sample and by gene from entire maf
  samples = factor(as.character(maf[,sample.column]), levels = sample.list)
  genes = factor(as.character(maf[,gene.column]), levels = gene.list)
  variants = as.character(maf[,variant.column])
 
  total.counts.by.sample = xtabs(~ variants + samples)
  total.counts.by.gene = xtabs(~ variants + genes)  

  ## Subset maf and get new info  
  maf = subset(maf, maf[,sample.column] %in% sample.list & maf[,gene.column] %in% gene.list)
  samples = factor(as.character(maf[,sample.column]), levels = sample.list)
  genes = factor(as.character(maf[,gene.column]), levels = gene.list)
  variants = as.character(maf[,variant.column])
  
  ## Perform category name substitution  							
  variant.id = rep(0, length(variants))
  for(i in 1:nrow(variant.categories)) {
    ind = grepl(variant.categories[i,"maf.classification"], variants)
    variant.id[ind] = variant.categories[i,"id"]
  }
  variant.id = as.integer(variant.id)
  
  ## If there are other categories in the maf that were not in the variant.categories variable, they will be listed
  variants.extra = unique(variants[variant.id == 0])
  if(length(variants.extra) > 0) {
    if(verbose == TRUE) { cat("The following variant categories were in the maf but not in the variant.categories variable:", "\n", paste(variants.extra, collapse="\n"), "\nThese will not be considered in the final matrix\n") }
  }
  
  variant.agg <- aggregate(formula = variant.id ~ genes + samples, FUN = max, na.rm=TRUE) 
  variant.matrix <- as.matrix(xtabs(formula = variant.id ~ genes + samples, data=variant.agg, drop.unused.levels = FALSE))
  class(variant.matrix) = "matrix"

  ## Fill in samples with missing data (i.e. NA samples)
  if(!is.null(na.samples)) {
    variant.matrix = cbind(variant.matrix, matrix(NA, nrow = nrow(variant.matrix), ncol = length(na.samples), dimnames = list(rownames(variant.matrix), na.samples)))
  }                      
  if(!is.null(na.genes)) {
    variant.matrix = rbind(variant.matrix, matrix(NA, nrow = length(na.genes), ncol = ncol(variant.matrix), dimnames = list(na.genes, colnames(variant.matrix))))
  }                      

  ## Order rows 
  if(nrow(variant.matrix) > 1) { 
    o = matrix.sort(variant.matrix, ordering=gene.order, margin="rows")
    variant.matrix = variant.matrix[o,,drop=FALSE]
  }  

  ## Order columns
  if(ncol(variant.matrix) > 1) {
    o = matrix.sort(variant.matrix, ordering=sample.order, margin="columns")
    variant.matrix = variant.matrix[,o,drop=FALSE]
  }

  return(list(matrix=variant.matrix, type="categorical", mapping=variant.categories, counts.by.gene=total.counts.by.gene, counts.by.sample=total.counts.by.sample))
}




is.overlapping = function(set1, set2) {
  ## Determines if two vectors contain the same values (not necessarily in the same order)
  flag = FALSE
  l1 = length(set1)
  l2 = length(set2)
  if(l1 == l2 & length(intersect(set1, set2)) == l1) {
    flag = TRUE
  }
  return(flag)
}


read.maf = function(filename, stringsAsFactors=FALSE) {
  maf = read.table(filename, header=TRUE, stringsAsFactors=stringsAsFactors, sep="\t", quote="")
  return(maf)
}


my.mapping = function() {
    require(RColorBrewer)
    variant.categories = data.frame(rbind(c(9, "Nonsense_Mutation", "Nonsense"),
								c(8, "Frame_Shift.*", "Frame Shift"),
								c(7, "[iI]n_[fF]rame.*", "In Frame Indel"),
								c(6, "Splice_Site", "Splice Site"),
								c(5, "Missense_Mutation", "Missense"),
								c(4, "Nonstop_.*|De_novo_Start.*|De_novo_Stop.*|Start_.*|Stop_.*|Translation_.*|Read\\-through.*", "Other"),
								c(3, "3'\\-?UTR|5'\\-?UTR|3'\\-?Flank|5'\\-?Flank|Intron", "UTR/Intron/Flanking"),
								c(2, "IGR", "IGR"),
								c(1, "RNA|lincRNA", "RNA")), stringsAsFactors=FALSE)
								
	return(variant.categories)	
}


tcga.mapping = function() {
  
  require(RColorBrewer)
  
  ## These 7 Categories are used in many of the TCGA papers
  id = c(-1, 1:6)
  maf.classification = c("Silent",
  						"Nonstop_.*|De_novo_Start.*|De_novo_Stop.*|Start_.*|Stop_.*|Translation_.*|Read\\-through.*",
  						"Missense_.*",
  						"Splice_Site.*",
  						"[iI]n_[fF]rame.*",
  						"Frame_Shift.*",
  						"Nonsense_.*")
  
  final.classification = c("Syn.", "Other non syn.", "Missense", "Splice site", "In frame indel", "Frame shift", "Nonsense")						
  
  brewer.pal.set1 = brewer.pal(7, "Set1")
  bg.col = c(brewer.pal(7, "Set1")[c(3, 1, 2, 4, 6, 5, 1)]) 
  
  return(data.frame(id, maf.classification, final.classification, bg.col, stringsAsFactors=FALSE))
}




tabler = function(tabler.obj, plot.samples="all", plot.samples.index=NULL,
  xlab = "", ylab = "", 
  horizontal.border = TRUE,
  vertical.border = TRUE,
  horizontal.border.args = list(lwd=0.25, col="grey90"),
  vertical.border.args = list(lwd=0.25, col="grey90"),
  box = TRUE,
  box.args = list(col="black", lwd=1),
  row.label = TRUE,
  column.label = TRUE,
  na.col = "grey90",
  bg.col = "white",
  row.freq = TRUE,
  row.label.args = list(),
  column.label.args = list(),
  row.freq.args = list(),
  image.args = list(),
  column.label.side = 'bottom', row.label.side = "left", row.freq.side = "right",
  y.adj = 0, # adjustment for row labels
  x.adj = 0, # adjustment for col labels
  ...)
  {

    mat = tabler.obj$matrix
    tabler.mapping = tabler.obj$mapping


    if(nrow(mat) > 1) {
      mat = mat[nrow(mat):1,]
    }
    row.freq.percent = round((rowSums(mat>0, na.rm=TRUE) / rowSums(!is.na(mat)))*100)
    row.freq.label = paste(row.freq.percent, "%", sep="")
    row.freq.label[row.freq.percent < 1] = "<1%"

    ## Remove unaltered samples if requested
    if(plot.samples == "mutated") {
      one.mut.present = colSums(mat == 0 | is.na(mat) | is.infinite(mat), na.rm=TRUE) != nrow(mat)
      mat = mat[,one.mut.present]
      print(plot.samples)
    } else if (plot.samples == "selected" & !is.null(plot.samples.index)) {
      mat = mat[,plot.samples.index]
    }


    tabler.mapping.id = as.integer(tabler.mapping$id)
    tabler.mapping.color = as.character(tabler.mapping$bg.col)
 
    ## Add color for zero if it isn't already in there
    if(sum(tabler.mapping.id == 0) == 0) {
      tabler.mapping.id = c(tabler.mapping.id, 0)
      tabler.mapping.color = c(tabler.mapping.color, bg.col)
    }
    ## Add color for NA samples if they exist in the matrix
    if(sum(is.na(mat)) > 0) {
      na.id = min(tabler.mapping.id)-1
     
      tabler.mapping.id = c(tabler.mapping.id, na.id)
      tabler.mapping.color = c(tabler.mapping.color, na.col)
      
      mat[is.na(mat)] = na.id
    }

    ## Order mapping ids and colors for use in image()
    ind = order(tabler.mapping.id)
    tabler.mapping.id = tabler.mapping.id[ind]
    tabler.mapping.color = tabler.mapping.color[ind]
    
    ## Remove values that were in the mapping scheme but not in the actual matrix
    id.in.matrix = tabler.mapping.id %in% c(mat)
    tabler.mapping.id = tabler.mapping.id[id.in.matrix]
    tabler.mapping.color = tabler.mapping.color[id.in.matrix]
    
    ## Rescale IDs such that values in matrix are increasing by one    
    ## Needed for image(), otherwise it behaves strange when trying to color
    new.mat = mat
    new.ids = 1:length(tabler.mapping.id)
    for(i in 1:length(new.ids)) {
      new.mat[mat == tabler.mapping.id[i]] = new.ids[i] 
    }
    
    ## Make matrix plot    
    image.required.args = list(1:ncol(new.mat), 1:nrow(new.mat), t(new.mat), col=tabler.mapping.color, axes=FALSE, xlab=xlab, ylab=ylab, zlim=c(min(new.mat, na.rm=TRUE)-1, max(new.mat, na.rm=TRUE)+1))
    image.args = c(image.required.args, image.args)
    do.call("image", image.args)
    
    ## Add Vertical borders
    if(vertical.border == TRUE) {
      vertical.segs = list(.5 + 1:(ncol(new.mat)-1), .5, .5 + 1:(ncol(new.mat)-1), .5 + nrow(new.mat))
      do.call("segments", c(vertical.segs, vertical.border.args))
    }
    
    ## Add Horizontal borders
    if(horizontal.border == TRUE) {
      horizontal.segs = list(.5, .5 + 1:(nrow(new.mat)-1), .5 + ncol(new.mat), .5 + 1:(nrow(new.mat)-1))
      do.call("segments", c(horizontal.segs, horizontal.border.args))
    }
    
    ## Add column labels
    if(column.label == TRUE) {
      if (column.label.side == "top") {
        x.side = 3
      } else {
        x.side = 1
      } 
      
      column.label.required.args = list(side=x.side, at = 1:ncol(new.mat), labels = colnames(new.mat), tick = FALSE, las = 2)
      do.call("axis", c(column.label.required.args, column.label.args))
    }  
    
    ## Add row labels
    if(row.label == TRUE) {
      if (row.label.side == "left") {
        y.side = 2
      } else {
        y.side = 4
      } 
      
      row.label.required.args = list(side=y.side, at = 1:nrow(new.mat), labels = rownames(mat), tick = FALSE, las = 1)
      do.call("axis", c(row.label.required.args, row.label.args))
    }
    
    ## Add row frequecies
    if(row.freq == TRUE) {
          
      if (row.freq.side == "left") {
        y.side = 2
      } else {
        y.side = 4
      } 
      row.freq.required.args = list(side=y.side, at = 1:nrow(new.mat), labels = row.freq.label, tick = FALSE, las = 1)
      do.call("axis", c(row.freq.required.args, row.freq.args)) 
    }  
    
    ## Plot box around matrix    
    if (box == TRUE) {
      do.call("box", c(list(which="plot"), box.args))
    }
   
  }







matrix.sort = function(M, ordering="frequency", margin="rows", binary.cutoff=0, decreasing=TRUE) {

  ## The rest of the function assumes ordering will be done by row. Transpose matrix if columns need sorting
  if(margin == "columns") {
    M = t(M)
  }
  
  ## Handle NAs as lowested possible kind of mutation
  M.min = min(M, na.rm=TRUE)
  if(M.min > -1) { M.min = -1 }
  M[is.na(M) | is.infinite(M)] = M.min

  ## Generate different matricies 
  index = 1:nrow(M)
  M.pos = M>binary.cutoff
  M.any = M != 0
  
  ## Order samples by overall frequency, mutation type, or a given set of indices or labels
  if(ordering == "frequency") {
    index = order(rowSums(M.pos, na.rm=TRUE), decreasing=decreasing)
  } else if (ordering == "mutation.type") {
    ## Sorts columns first by whether or not they are mutated and then by the type of mutation
    M.combined = cbind(M.pos, M.any, M)
    index = do.call("order", c(lapply(1:ncol(M.combined), function(x) M.combined[,x]), decreasing=decreasing))
  } else if (class(ordering) == "character" & is.overlapping(ordering, rownames(M))) {
    index = match(rownames(M), ordering) 
  } else if (class(ordering) == "integer" & is.overlapping(ordering, 1:nrow(M))) {
    index = ordering
  } else if (sample.order != "none") {
    warning('ordering needs to be one of "none", "frequency", "mutation.type", a character vector matching the names of the matrix, or an integer vector matching the number of elements of the matrix margin. Ordering of samples not performed')
  }

  return(index)
}






multi.tabler = function(table.list, sample.order="mutation.type", close.screens=TRUE, between.table.space = 0.025, row.space = 0.025, mar=c(4.1,4.1), bg.col="white", na.col="grey90", labCol=TRUE, plot.samples="all", ...) {
  
  ntable = length(table.list)
  
  ## Check to verify that column names match for each table in the list  
  table.colnames = colnames(table.list[[1]][["matrix"]])
  if(ntable > 1) {
    for(i in 2:ntable) {
      temp.table = table.list[[i]][["matrix"]]
      if(!is.overlapping(table.colnames, colnames(temp.table))) {
        stop(paste("Column names in table", i, "do not match those from table 1. Cannot align columns tables properly\n", sep=" "))
      }
    }
  }
  
  ## Create meta-matrix for ordering, start with last one
  meta.matrix = table.list[[ntable]][["matrix"]]
  meta.matrix = meta.matrix[,table.colnames,drop=FALSE]

  
  if(ntable > 1) {
    for(i in (ntable-1):1) {

      ## From the last table to the first, increase or decrease the ids so they are non-overlapping 
      ## with all other tables. Then concatenate into one super table for sorting
      meta.matrix.max = max(meta.matrix, na.rm=TRUE)
      meta.matrix.min = min(meta.matrix, na.rm=TRUE)
      matrix.to.add = table.list[[i]][["matrix"]][,table.colnames,drop=FALSE]
      

      ## Shift values to be greater/less than the values currently in the meta.matrix
      ## This will give priority to tables earlier in the table.list for sorting
      ind = matrix.to.add < 0 & !is.na(matrix.to.add)
      matrix.to.add[ind] = matrix.to.add[ind] + meta.matrix.min
      
      ind = matrix.to.add > 0 & !is.na(matrix.to.add)
      matrix.to.add[ind] = matrix.to.add[ind] + meta.matrix.max
      
      meta.matrix = rbind(matrix.to.add[,table.colnames,drop=FALSE], meta.matrix)
    }
  } 

  ## Determine which samples do and do not have at least one mutation
  if(plot.samples == "all") {
    samples.to.show = rep(TRUE, ncol(meta.matrix))
  } else {
    samples.to.show = rep(FALSE, ncol(meta.matrix))
    for(i in 1:ntable) {
      is.mut = apply(table.list[[i]][["matrix"]][,table.colnames, drop=FALSE] != 0, 2, sum, na.rm=TRUE) > 0
      samples.to.show[is.mut] = TRUE  
    }
  } 

  ## Perform ordering of samples across all tablers
#  meta.matrix.to.show = meta.matrix[,samples.to.show]
#  global.sample.order = matrix.sort(meta.matrix.to.show, margin="columns", ordering=sample.order)
  global.sample.order = matrix.sort(meta.matrix, margin="columns", ordering=sample.order)
  global.sample.names = colnames(meta.matrix)[global.sample.order]

  samples.to.show.order = samples.to.show[global.sample.order]
  
  
  ## Plot tables one by one
  table.nrows = unlist(lapply(1:ntable, function(i) { nrow(table.list[[i]][["matrix"]]) }))
  total.nrow = sum(table.nrows)
  
  total.table.space = row.space*table.nrows
  between.table.shift = c(cumsum(rep(between.table.space, ntable)), 0) ## adjustment for space between tables
  
  ss.top = c(1, 1-cumsum(total.table.space)) - between.table.shift 
  ss.bottom = bottom = c(1-cumsum(total.table.space), 0) - between.table.shift
  ss.left = 0  ## par(mar) will be used for now to move the sides of the plots in
  ss.right = 1

  ss.matrix = cbind(ss.left, ss.right, ss.bottom, ss.top)
  
  ## Remove dummy row at the end and ensure it is a matrix (for the case when only one element in table.list)
  ss.matrix = matrix(ss.matrix[-nrow(ss.matrix),], byrow=FALSE, ncol=4) 

  ## Minimum amount of space at the bottom for table column labels
  min.bottom.space = 0.1
  if(min(ss.matrix[,3]) < min.bottom.space) {
    ss.matrix[,3:4] = ss.matrix[,3:4] + min.bottom.space + abs(min(ss.matrix[,3]))
  }
  
  
  ## Check to see if ss.matrix values are outside [0,1] interval and adjust if necessary
  ss.matrix = rescale.screen(ss.matrix)
  
  screen.index = split.screen(ss.matrix)

  plot.col.labels = FALSE
  for(i in 1:ntable) {
    
    if(i == ntable & labCol == TRUE) {
      plot.col.labels = TRUE
    }
    
    screen(screen.index[i])
    tabler.matrix = table.list[[i]][["matrix"]][,global.sample.names,drop=FALSE]
    tabler.mapping = table.list[[i]][["mapping"]]
    new.tabler = list(matrix=tabler.matrix, mapping=tabler.mapping)
    
    ## Plot table
    par(mar=c(0,mar[1],0,mar[2]), mgp=c(0,0.1,0.1))
    tabler(new.tabler, column.label=plot.col.labels, plot.samples="selected", plot.samples.index=samples.to.show.order, bg.col=bg.col, na.col=na.col, ...)

  }  

  if(close.screens == TRUE) {
    close.screen(screen.index)
  }
  return(meta.matrix)
}



rescale.screen = function(screen.matrix) {
  if(ncol(screen.matrix) != 4) {
    stop("Matrix for split.screen must be four columns")
  }

  ## Rescale left and right boundries if they excede [0,1]
  if(min(screen.matrix[,1]) < 0) {
    screen.matrix[,1:2] = screen.matrix[,1:2] + abs(min(screen.matrix[,1]))
  }
  if(max(screen.matrix[,4]) > 1) {
    screen.matrix[,1:2] = screen.matrix[,1:2] / max(screen.matrix[,2])
  }
  
  ## Rescale top and bottom boundries if they excede [0,1]
  if(min(screen.matrix[,3]) < 0) {
    screen.matrix[,3:4] = screen.matrix[,3:4] + abs(min(screen.matrix[,3]))
  }
  if(max(screen.matrix[,4]) > 1) {
    screen.matrix[,3:4] = screen.matrix[,3:4] / max(screen.matrix[,4])
  }

  return(screen.matrix)
}
