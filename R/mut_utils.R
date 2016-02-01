

## Examples
## "c.326+2_326+3delTA"
## "c.1316A>C"
## c.1509_1510insGCCTAT

parse_CDS_str = function(str) {

  str.regex = "^c.([0-9(]+)([+|-][1-9?]*){0,1}_{0,1}([0-9)]*)([+|-][1-9?]*){0,1}([ACGTN]*)(>|del|ins)([ACGTN]*)"
  cds.match = str_match(str, str.regex)
  
  no.match = sum(is.na(cds.match[,1]))
  if(no.match > 0) {
    no.match.str = paste(head(unique(str[is.na(cds.match[,1])]), 5), collapse=",")
    warn.str = paste(no.match, " CDS strings not properly matched and parsed. Example(s): ", no.match.str, sep="")
    warning(warn.str)
  }
  
  colnames(cds.match) = c("CDS", "Exonic_start", "Relative_exon_start", "Exonic_end", "Relative_exon_end", "Reference_Allele", "Type", "Alternate_Allele")
  return(cds.match)
}  








## Function to test each row of a matrix with the FET
fet = function(m, label, alternative = c("two.sided", "greater", "less"), reorder=TRUE) {
  require(reshape2)
  
  alternative = match.arg(alternative)
  if(!is.factor(label)) { label = as.factor(label) }

  pval = c()
  OR = c()
  freq = c()
  tot = c()
  for(i in 1:nrow(m)) {
    mut = factor(ifelse(m[i,] > 0, "Mutated", "Wt"), levels=c("Wt", "Mutated"))

    f = fisher.test(label, mut, alternative=alternative)
    pval = c(pval, f$p.value)
    OR = c(OR, ifelse(!is.null(f$estimate), f$estimate, NA))

    ## Summary stats for each group
	ta = table(mut, label)
	ta.melt = melt(ta)
    ta.melt.final = ta.melt[,3]
    names(ta.melt.final) = paste("Total", ta.melt[,1], ta.melt[,2], sep="_")

	ta.frac = sweep(ta, 2, apply(ta, 2, sum), "/")
    ta.frac.melt = melt(ta.frac)
    ta.frac.melt.final = ta.frac.melt[,3]
    names(ta.frac.melt.final) = paste("Fraction", ta.frac.melt[,1], ta.frac.melt[,2], sep="_")
    
    freq = rbind(freq, ta.frac.melt.final)
    tot = rbind(tot, ta.melt.final)    
  }
    

  fdr = p.adjust(pval, "fdr")
  
  
  final = data.frame(Gene=rownames(m), tot, freq, OR=OR, Pvalue=pval, FDR=fdr, row.names=rownames(m), check.names=FALSE, stringsAsFactors=FALSE)
  
  if(reorder==TRUE) {
    o = order(pval)   
    final = final[o,]
  }
  
  return(final)  
}



