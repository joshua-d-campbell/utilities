#################
# read.maf
#
# Read in a MAF file into a data.frame, data.table, or GenomicRanges object
#
#################
read.maf = function(file, class=c("dt", "df", "gr"), keys=c("Chromosome", "Start_position", "End_position", "Strand", "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2"), standardize=TRUE, showProgress=FALSE, ...) {
	require(data.table)
	require(GenomicRanges)
    
	class = match.arg(class)
	
	if(class == "df") {
	
	  temp.maf = suppressWarnings(fread(file, data.table=FALSE, showProgress=showProgress, ...))
  	  if(standardize == TRUE) { standardize.maf.columns(temp.maf) }
  	  
	} else if (class == "dt" | class == "gr") {
	
	  temp.maf = suppressWarnings(fread(file, data.table=TRUE, showProgress=showProgress, ...))
  	  if(standardize == TRUE) { standardize.maf.columns(temp.maf) } ## Standardize before setting keys
	  setkeyv(temp.maf, keys)	  
	  
    }
    
    if(class == "gr") {
      maf.gr = GRanges(seqnames=Rle(temp.maf$Chromosome), IRanges(names=temp.maf$Chromosome, start=temp.maf$Start_position, end=temp.maf$End_position), strand=temp.maf$Strand, temp.maf)
      return(maf.gr)
    } else {
      return(temp.maf)
    }
}



##################
# standarize.maf.columns
#
# Some programs such as MutSig will change the column names of the maf. This function automatically checks and fixes the columns
#
##################
standardize.maf.columns = function(maf) {
  colnames.fix = c("Start_Position", "End_Position")
  colnames.correct = c("Start_position", "End_position")

  colnames.match = match(colnames.fix, colnames(maf))
  ix = !is.na(colnames.match)
  if(sum(ix) > 0) { colnames(maf)[colnames.match[ix]] = colnames.correct[ix] }
  return(maf)
}




##################
# maf.coding
#
# Scans "Variant_Classification" field in maf and outputs a subseted maf or an index correspoding to the
# rows in the maf that have a coding mutation
#
##################
maf.coding = function(maf, index=FALSE)
{  
  ix = maf$Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Stop_Codon_Del", "Stop_Codon_Ins", "Start_Codon_Del", "Start_Codon_Ins")	

  if(index == TRUE) {
    return(ix)	  
  } else {
    return(maf[ix,])
  }
}



#################
# maf.exonic
#
# Scans "Variant_Classification" field in maf and outputs a subseted maf or an index correspoding to the
# rows in the maf that have a coding or a silent mutation
#
##################
maf.exonic = function(maf, index=FALSE)
{  
  ix = maf$Variant_Classification %in% c("Silent", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Stop_Codon_Del", "Stop_Codon_Ins", "Start_Codon_Del", "Start_Codon_Ins")	

  if(index == TRUE) {
    return(ix)	  
  } else {
    return(maf[ix,])
  }
}



#################
# maf.genic
#
# Scans "Variant_Classification" field in maf and outputs a subseted maf or an index correspoding to the
# rows in the maf that have a coding, silent, UTR, or intronic mutation
#
##################
maf.genic = function(maf, index=FALSE)
{  
  ix = maf$Variant_Classification %in% c("3'UTR", "5'UTR", "Intron", "Silent", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site", "Splice_Site", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Stop_Codon_Del", "Stop_Codon_Ins", "Start_Codon_Del", "Start_Codon_Ins")	

  if(index == TRUE) {
    return(ix)	  
  } else {
    return(maf[ix,])
  }
}


###################
# mutation.context.snv192
#
# Summarizes the 192 mutation contexts per sample
#
##################

mutation.context.snv192 = function(maf) {

  require(Biostrings)
  
  final.mut.type = rep(NA, nrow(maf))
  final.mut.context = rep(NA, nrow(maf))
  final.mut.strand = rep(NA, nrow(maf))
    
  ## Get mutation type
  initial.maf.type = paste(maf$Reference_Allele, ">", maf$Tumor_Seq_Allele2, sep="")
  
  ## Get mutation context info for those on "+" strand
  forward.change = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  ind = maf$Variant_Type == "SNP" & initial.maf.type %in% forward.change

  final.mut.type[ind] = paste(maf$Reference_Allele[ind], ">", maf$Tumor_Seq_Allele2[ind], sep="")
  final.mut.context[ind] = toupper(substring(as.character(maf$ref_context[ind]), 10,12))
  final.mut.strand[ind] = ifelse(maf$Transcript_Strand[ind] == "+", "T", "U")

  ## Get mutation context info for those on "-" strand
  rev.change = c("A>G","A>T","A>C","G>T","G>C","G>A")
  ind = maf$Variant_Type == "SNP" & initial.maf.type %in% rev.change

  ## Reverse complement the context so only 6 mutation categories instead of 12
  rev.context = reverseComplement(DNAStringSet(maf$ref_context[ind]))
  rev.refbase = reverseComplement(DNAStringSet(maf$Reference_Allele[ind]))
  rev.altbase = reverseComplement(DNAStringSet(maf$Tumor_Seq_Allele2[ind]))

  final.mut.type[ind] = paste(as.character(rev.refbase), ">", as.character(rev.altbase), sep="")
  final.mut.context[ind] = toupper(substring(as.character(rev.context), 10,12))
  final.mut.strand[ind] = ifelse(maf$Transcript_Strand[ind] == "-", "T", "U")
  
  
  maf.mut.id = paste(final.mut.type, final.mut.context, final.mut.strand, sep="_")
  Tumor_ID = as.factor(maf$Tumor_Sample_Barcode)

  ## Define all mutation types for 196 substitution scheme
  b1 = rep(rep(c("A", "C", "G", "T"), each=24), 2)
  b2 = rep(rep(c("C", "T"), each=12), 8)
  b3 = rep(c("A", "C", "G", "T"), 48)
  mut.trinuc = apply(cbind(b1, b2, b3), 1, paste, collapse="")
  mut.type = rep(rep(rep(forward.change, each=4), 4), 2)
  mut.strand = rep(c("T", "U"), each=96)

  mut.id = apply(cbind(mut.type, mut.trinuc, mut.strand), 1, paste, collapse="_")  

  Mutation = factor(maf.mut.id, levels=mut.id)
  
  mut.table = xtabs(~ Mutation + Tumor_ID)
  mut.summary = data.frame(Mutation, Type=final.mut.type, Context=final.mut.context, Strand=final.mut.strand, stringsAsFactors=FALSE)
  
  return(list(mutation_table=mut.table, maf_mutations=mut.summary))  
}
