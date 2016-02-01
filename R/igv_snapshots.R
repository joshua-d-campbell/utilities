maf2igvcommand = function(pair.info, sample.info, maf, outfile, bams=c("clean_bam_file_capture", "clean_bam_file_wgs", "bam_file_rna_analysis_ready"), snapshot.dir, snapshot.name=c("Tumor_Sample_Barcode", "Normal_Sample_Barcode", "Chromosome", "Start_Position", "End_position", "Hugo_Symbol", "Variant_Classification", "Protein_Change", "Reference_Allele", "Tumor_Seq_Allele2"), genome="hg19", max.height=400, window=20) {
  unlink(outfile)
  con = file(outfile, "w")
    
  ## Check for bam columns in sample.info
  bad.columns = setdiff(bams, colnames(sample.info))
  if(length(bad.columns) > 0) {
    cat("The following bam columns were not found in the sample file and will not be given snapshots:", paste(bad.columns, collapse=", "), "\n")
  }
  bams = intersect(bams, colnames(sample.info))

  bad.columns = setdiff(snapshot.name, colnames(maf))
  if(length(bad.columns) > 0) {
    cat("The following columns were not found in the sample file and will not be used in the shapshot name:", paste(bad.columns, collapse=", "), "\n")
  }
  snapshot.name = intersect(snapshot.name, colnames(maf))
  
  ## Convert any factors in maf to character
  i <- sapply(maf, is.factor)
  maf[i] <- lapply(maf[i], as.character)
  
  ## Interate through each pair and set up snapshot text
  for(i in 1:nrow(pair.info)) {
  
    pair.id = pair.info[i,"pair_id"]
    tumor.id = pair.info[i,"case_sample"]
    normal.id = pair.info[i, "control_sample"]
    
    tumor.sample.info = subset(sample.info, sample_id == tumor.id)
    normal.sample.info = subset(sample.info, sample_id == normal.id)

    maf.s = subset(maf, Tumor_Sample_Barcode == tumor.id & Matched_Norm_Sample_Barcode == normal.id)
    
    if(nrow(tumor.sample.info) == 0) {
      cat("No tumor sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(tumor.sample.info) > 1) {
      cat("More than one tumor sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(normal.sample.info) == 0) {
      cat("No normal sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(normal.sample.info) > 1) {
      cat("More than one normal sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(maf.s) == 0) {
      cat("No mutations found for pair", pair.id, "... Skipping", "\n")
      next()
    } else {
      cat("Generating snapshots for", nrow(maf.s), "mutations for pair", pair.id, "\n")
    }
    
    ## Set up global variables for snapshots
    cat("new", "\n", file=con)
    cat("genome", genome, "\n", file=con)
    cat("maxPanelHeight", max.height, "\n", file=con)
    cat("snapshotDirectory", snapshot.dir, "\n", file=con)
    
    ## Print text for loading bam files
    
    for(j in bams) {
      tumor.bam = tumor.sample.info[1,j]
      normal.bam = normal.sample.info[1,j]
      
      if(tumor.bam != "" & !is.na(tumor.bam)) {
        cat("load", tumor.bam, "\n", file=con)
      } else {
        cat(j, "not found for", tumor.id, "\n")
      }
      if(normal.bam != "" & !is.na(normal.bam)) {
        cat("load", normal.bam, "\n", file=con)
      } else {
        cat(j, "not found for", normal.id, "\n")
      }
    }
    cat("echo loaded", "\n", file=con)
    
    ## Print text for snapshots of each position
    for(j in 1:nrow(maf.s)) {
      cat("goto", paste(maf.s$Chromosome[j], ":", maf.s$Start_position[j]-window, "-", maf.s$Start_position[j]+window, sep=""), "\n", file=con)
      cat("sort base", "\n", file=con)
      
      png.filename = paste(paste(maf.s[j,snapshot.name], collapse="_"), ".png", sep="")
      cat("snapshot", png.filename, "\n", file=con)
    }
  }
  
  cat("exit", "\n", file=con)
}








positions2igvcommand = function(pair.info, sample.info, positions, outfile, bams=c("clean_bam_file_capture", "clean_bam_file_wgs", "bam_file_rna_analysis_ready"), snapshot.dir, genome="hg19", max.height=400, window=30) {
  unlink(outfile)
  con = file(outfile, "w")
  
  ## Check for bam columns in sample.info
  bad.columns = setdiff(bams, colnames(sample.info))
  if(length(bad.columns) > 0) {
    cat("The following bam columns were not found in the sample file and will not be given snapshots:", paste(bad.columns, collapse=", "), "\n")
  }
  bams = intersect(bams, colnames(sample.info))
  
  ## Interate through each pair and set up snapshot text
  for(i in 1:nrow(pair.info)) {
  
    pair.id = pair.info[i,"pair_id"]
    tumor.id = pair.info[i,"case_sample"]
    normal.id = pair.info[i, "control_sample"]
    
    tumor.sample.info = subset(sample.info, sample_id == tumor.id)
    normal.sample.info = subset(sample.info, sample_id == normal.id)

    if(nrow(tumor.sample.info) == 0) {
      cat("No tumor sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(tumor.sample.info) > 1) {
      cat("More than one tumor sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(normal.sample.info) == 0) {
      cat("No normal sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(normal.sample.info) > 1) {
      cat("More than one normal sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else {
      cat("Generating snapshots for", nrow(positions), "mutations for pair", pair.id, "\n")
    }
    
    ## Set up global variables for snapshots
    cat("new", "\n", file=con)
    cat("genome", genome, "\n", file=con)
    cat("maxPanelHeight", max.height, "\n", file=con)
    cat("snapshotDirectory", snapshot.dir, "\n", file=con)
    
    ## Print text for loading bam files
    
    for(j in bams) {
      tumor.bam = tumor.sample.info[1,j]
      normal.bam = normal.sample.info[1,j]
      
      if(tumor.bam != "" & !is.na(tumor.bam)) {
        cat("load", tumor.bam, "\n", file=con)
      } else {
        cat(j, "not found for", tumor.id, "\n")
      }
      if(normal.bam != "" & !is.na(normal.bam)) {
        cat("load", normal.bam, "\n", file=con)
      } else {
        cat(j, "not found for", normal.id, "\n")
      }
    }
    cat("echo loaded", "\n", file=con)
    
    ## Print text for snapshots of each position
    for(j in 1:nrow(positions)) {
      cat("goto", paste(positions$Chromosome[j], ":", positions$Start_Position[j]-window, "-", positions$Start_Position[j]+window, sep=""), "\n", file=outfile, append=TRUE)
      cat("sort base", "\n", file=con)
      
      png.filename = paste(paste(c(tumor.id, normal.id, positions[j,]), collapse="_"), ".png", sep="")
      cat("snapshot", png.filename, "\n", file=con)
    }
  }
  
  cat("exit", "\n", file=con)
  close(con)
}







maf_multi2igvcommand = function(sample.mapping, maf, outfile, bams=c("clean_bam_file_capture", "clean_bam_file_wgs", "bam_file_rna_analysis_ready"), snapshot.dir, snapshot.name=c("Tumor_Sample_Barcode", "Normal_Sample_Barcode", "Chromosome", "Start_Position", "End_position", "Hugo_Symbol", "Variant_Classification", "Protein_Change", "Reference_Allele", "Tumor_Seq_Allele2"), genome="hg19", max.height=400, window=20) {
  unlink(outfile)
  con = file(outfile, "w")
    
  ## Check for bam columns in sample.info
  bad.columns = setdiff(bams, colnames(sample.info))
  if(length(bad.columns) > 0) {
    cat("The following bam columns were not found in the sample file and will not be given snapshots:", paste(bad.columns, collapse=", "), "\n")
  }
  bams = intersect(bams, colnames(sample.info))

  bad.columns = setdiff(snapshot.name, colnames(maf))
  if(length(bad.columns) > 0) {
    cat("The following columns were not found in the sample file and will not be used in the shapshot name:", paste(bad.columns, collapse=", "), "\n")
  }
  snapshot.name = intersect(snapshot.name, colnames(maf))
  
  ## Convert any factors in maf to character
  i <- sapply(maf, is.factor)
  maf[i] <- lapply(maf[i], as.character)
  
  ## Interate through each pair and set up snapshot text
  for(i in 1:nrow(pair.info)) {
  
    pair.id = pair.info[i,"pair_id"]
    tumor.id = pair.info[i,"case_sample"]
    normal.id = pair.info[i, "control_sample"]
    
    tumor.sample.info = subset(sample.info, sample_id == tumor.id)
    normal.sample.info = subset(sample.info, sample_id == normal.id)

    maf.s = subset(maf, Tumor_Sample_Barcode == tumor.id & Normal_Sample_Barcode == normal.id)
    
    if(nrow(tumor.sample.info) == 0) {
      cat("No tumor sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(tumor.sample.info) > 1) {
      cat("More than one tumor sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(normal.sample.info) == 0) {
      cat("No normal sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(normal.sample.info) > 1) {
      cat("More than one normal sample found for pair", pair.id, "... Skipping", "\n")
      next()
    } else if (nrow(maf.s) == 0) {
      cat("No mutations found for pair", pair.id, "... Skipping", "\n")
      next()
    } else {
      cat("Generating snapshots for", nrow(maf.s), "mutations for pair", pair.id, "\n")
    }
    
    ## Set up global variables for snapshots
    cat("new", "\n", file=con)
    cat("genome", genome, "\n", file=con)
    cat("maxPanelHeight", max.height, "\n", file=con)
    cat("snapshotDirectory", snapshot.dir, "\n", file=con)
    
    ## Print text for loading bam files
    
    for(j in bams) {
      tumor.bam = tumor.sample.info[1,j]
      normal.bam = normal.sample.info[1,j]
      
      if(tumor.bam != "" & !is.na(tumor.bam)) {
        cat("load", tumor.bam, "\n", file=con)
      } else {
        cat(j, "not found for", tumor.id, "\n")
      }
      if(normal.bam != "" & !is.na(normal.bam)) {
        cat("load", normal.bam, "\n", file=con)
      } else {
        cat(j, "not found for", normal.id, "\n")
      }
    }
    cat("echo loaded", "\n", file=con)
    
    ## Print text for snapshots of each position
    for(j in 1:nrow(maf.s)) {
      cat("goto", paste(maf.s$Chromosome[j], ":", maf.s$Start_Position[j]-window, "-", maf.s$Start_Position[j]+window, sep=""), "\n", file=outfile, append=TRUE)
      cat("sort base", "\n", file=con)
      
      png.filename = paste(paste(maf.s[j,snapshot.name], collapse="_"), ".png", sep="")
      cat("snapshot", png.filename, "\n", file=con)
    }
  }
  
  cat("exit", "\n", file=con)
}








