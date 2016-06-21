## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
#
# General utility functions
#

MARCIN.COMPUTER = '7cc-7141-7010a-v100';
MARCIN.COMPUTER = '69.173.118.211';
MARCIN.COMPUTER = '69.173.118.172';
MARCIN.COMPUTER = 'mbr19d-4a1.local'
MARCIN.COMPUTER = 'vpn5-53'
MARCIN.COMPUTER = 'wm19d-4a1'
DEFAULT.HTAB.FILE = '~/public_html/htab.html';


# used for genome plot
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        {
                          if (class(x) == 'Hits')
                            c(queryLength(x), subjectLength(x))
                          else
                            as.numeric(dim(x))[1:2]
                        }))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

## shorthand listing largest objects in the workspace
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

## lapply(list, dim) shortcut
ldim = function(l)
  {
    return(lapply(l, function(x) {if (!is.null(dim(x))) dim(x) else length(x)}))
  }


############################
# qq_chisq 
#
# plots qq plot for observed chi squared vector obs against expected with df = df
#############################
qq_chisq = function(obs, df=1, highlight = c(), hexbins = NULL, samp = NULL)
  {
    if (!is.null(hexbins))
      require(hexbin)
    
    obs = obs[!is.na(obs)]

    if (!is.null(samp))
      if (samp<length(obs))
        obs = sample(obs, samp)
    
    exp = rchisq(length(obs), df = df)
    ord = order(obs)
    smax = max(c(exp, obs))
    colors = vector(mode = "character", length = length(obs)); colors[] = "black"; colors[highlight] = "red";
    colors = factor(colors[ord]);

    if (!is.null(hexbins))
      plot(hexbin(sort(exp), obs[ord], xbins = hexbins,  xlab = "Expected Chi^2", ylab = "Observed Chi^2"))
    else
      {
        plot(sort(exp), obs[ord], xlab = "Expected Chi^2", ylab = "Observed Chi^2", xlim = c(0, smax), ylim = c(0, smax), col = colors, pch = 18);
        lines(x=c(0, smax), y = c(0, smax), col = "gray");
      }
  }

#########
# file.name
#
# grabs filenames from list of paths
########
file.name = function(paths)
  {
    return(gsub('(^|(.*\\/))?([^\\/]*)', '\\3', paths))
  }

#########
# file.dir
#
# grabs file.dirs from list of paths
########
file.dir = function(paths)
  {
    return(gsub('(^|(.*\\/))?([^\\/]*)$', '\\2', paths))
  }

########
# splits a single string according to fixed widths contained in fw (ie each components i of fw denotes the width of field i in string str
########
strsplit.fwf = function(str, fw)
  {
    if (length(str)>1)
      {
        warning('String should be of length 1, only taking first element')
        str = str[1];
      }
    
    cs = cumsum(fw);
    return(substr(rep(str, length(fw)), cs-fw+1,c(cs[1:(length(cs)-1)], nchar(str))))
  }

########
# Returns vector of line counts for each file in path 
########
line.counts = function(paths)
  {
    out = rep(NA, length(paths))
    ix = which(file.exists(paths))
    out[ix] = sapply(paths, function(x) { p = pipe(paste('cat ', x, ' | wc -l ')); as.numeric(readLines(p)); close(p)});
    return(out)
  }


##########
# pk
#
# "pk"'s at vector or data frame ie samples n rows from it with replacement (or max rows / items if less than n)#
##########
pk = function(x, n = 10)
  {
    if (inherits(x, 'data.frame'))
      {
        n = min(nrow(x), n)
        return(x[sample(1:nrow(x), n), ])        
      }
    else
      {
        n = min(length(x), n)
        return(sample(x, n))
      }      
  }


##############
# border - orders rows of a logical / binary matrix treating each row as binary number with digits encoded as TRUE / FALSE values of entries
#
##############
border = function(B, na.rm = TRUE)
  {
    B = array(as.logical(B), dim = dim(B))
    tmp = vector(mode = "numeric", length = nrow(B));
    if (na.rm)
      B[is.na(B)] = FALSE;
    for (i in 1:ncol(B))
        tmp = tmp + 2^(ncol(B)-i)*as.numeric(B[,i]==1);
    return(order(tmp))
  }

#########
# Returns vector of column names for list of files, using sep as delimiter
#########
column.names = function(paths, sep = "\t")
  {
    out = list();
    out[paths] = NA;
    ix = file.exists(paths);
    out[ix] = strsplit(sapply(paths[ix], function(x) readLines(x, 1)), sep)
    return(out)
  }

############################
# qq_pval
#
# plots qq plot for observed pval vs uniform
#############################
qq_pval = function(obs, highlight = c(), title="", hexbins = NULL, samp = NULL, lwd = 1, bestfit=T, color='black', input.pch=18, input.cex=1, conf.lines=T, input.MAX=NULL, qvalues=NULL, genes=NULL, ...)
{	
	obs = -log10(obs[!is.na(obs)])
	obs = obs[!is.infinite(obs)]
	
	if (!is.null(samp))
		if (samp<length(obs))
			obs = sample(obs, samp)
	
	N <- length(obs)
	## create the null distribution 
	## (-log10 of the uniform)
	exp <- -log(1:N/N,10)
	
	if (is.null(input.MAX))
		MAX <- max(obs,exp) + 0.5
	else
		MAX <- input.MAX
	
	c95 <- rep(0,N)
	c05 <- rep(0,N)
	
	for(i in 1:N){
		c95[i] <- qbeta(0.975,i,N-i+1)
		c05[i] <- qbeta(0.025,i,N-i+1)
	}
	
	if (conf.lines){
		## plot the two confidence lines
		plot(exp, -log(c95,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", axes=FALSE, xlab="", ylab="")
		par(new=T)
		plot(exp, -log(c05,10), ylim=c(0,MAX), xlim=c(0,MAX), type="l", axes=FALSE, xlab="", ylab="")
		par(new=T)
		p1 <- c(exp[1], exp[1])
		p2 <- c(-log(c95,10)[1], -log(c05,10)[1])
		lines(x=p1, y=p2)
		x.coords <- c(exp,rev(exp))
		y.coords <- c(-log(c95,10),rev(-log(c05,10)))
		polygon(x.coords, y.coords, col='light gray', border=NA)
		par(new=T)
	}
	
	ord = order(obs)
	
	colors = vector(mode = "character", length = length(obs)); colors[] = "black"; colors[highlight] = "red";
	colors = factor(colors[ord]);
	
        dat = data.frame(x = sort(exp), y = obs[ord]);
        plot(dat$x, dat$y, xlab = expression(Expected -log[10](italic(P))), ylab = expression(Observed -log[10](italic(P))), xlim = c(0, MAX), ylim = c(0, MAX), pch=input.pch, cex=input.cex, bg=color, main=title, ...);
        if (!is.null(qvalues)){
          genes.to.label <- which(qvalues[ord] <= 0.2);
          genes <- genes[ord];
          if (length(genes.to.label > 0)){
            text(dat$x[genes.to.label], dat$y[genes.to.label], labels=genes[genes.to.label], pos=3);
          }
        }
        lines(x=c(0, MAX), y = c(0, MAX), col = "black", lwd = lwd);
        lambda = lm(y ~ x-1, dat)$coefficients;
        if (bestfit)
          {
            lines(x=c(0, MAX), y = c(0, lambda*MAX), col = "red", lty = 2, lwd = lwd);
            legend('bottomright',sprintf('lambda = %.2f', lambda), text.col='red', bty='n')
          }        
      }


##############################
# wfplot
#
# Quick waterfall plot
#
# data is a numeric vector
# labels are text labels of the same length as data
# col is either (1) an unamed list of unique colors (2) a named list mapping unique labels to colors
##############################
wfplot = function(data, labels = NULL, col = NULL, leg.pos = NULL, ...)
  {
    ix = order(data);
    ulab = unique(labels)

    if (!is.null(col))
      {
        if (is.null(names(col)))
          {
            if (length(ulab)>2)
              col = brewer.pal(length(ulab), 'Set3')
            else
              col = c('gray', 'red')

            col = col[match(labels, ulab)]
            names(col) = labels;
          }
        else
          {
            og.col = col;
            col = col[labels];
            ulab = intersect(names(og.col), names(col))
          }
      }
    
    barplot(data[ix], col = col[ix], border = FALSE, ...)

    if (is.null(leg.pos))
      leg.pos = c(mean(par('usr')[1:2])/4, 3*data[which.max(abs(data))]/4)
        
    legend(leg.pos[1], leg.pos[2], legend = ulab, fill = col[ulab])
  }
               
               


############
# Takes a character vector and makes an expression for creating that character vector
############
list.expr = function(x)
  {
    paste("c('", paste(x, sep = "", collapse = "', '"), "')", sep = "")
  }

########
# toggles options error recover / NULL
########
fuckr = function()
  {
    if (!is.null(options()$error))
      {
        options(error = NULL);
        print('Options error set to NULL');
      }
    else
      {
        options(error = recover);
        print('Options error set to recover');
      }
  }

#################
## flatten
##
## flattens 3rd dim of 3D array along cdim 1 (ie rows) or cdim 2 (ie cols) pasting together the appropriate combinations of dimnames with sep "sep"
## or if sep = NULL, then just dropping the 3rd dimension names
##
#################
flatten = function(A, cdim = 2, sep = "_")
{
  if (!(cdim==1 | cdim ==2))
    stop('cdim must be 1 or 2')

  ind = order(rep(c(1:dim(A)[cdim]), dim(A)[3]));
  
  out = A[,,1];

  if (cdim == 2)
    {
      if (dim(A)[3]>1)
        for (i in 2:dim(A)[3])
          out = cbind(out, A[,,i]);            
      dimnames(A)[[1]] = dimnames(A)[[1]]
    }
      
  if (cdim == 1)
    {
      if (dim(A)[3]>1)
        for (i in 2:dim(A)[3])
          out = rbind(out, A[,,i]);            
      dimnames(A)[[2]] = dimnames(A)[[2]]
    }

  out = out[,ind]; #reshuffle to get desired ordering  
  newdimnames = rep(dimnames(A)[[cdim]], each = dim(A)[3]);
  if (!is.null(sep))
    newdimnames = paste(newdimnames, dimnames(A)[[3]], sep = sep);  
  dimnames(out)[[cdim]] = newdimnames;

  return(out)
}

##################
# Makes bsub command that wraps shell command "cmd" to send to queue "queue"
# redirebmccting output / error etc streams to path prefixed by "jname",
# optional_args: maximum memory requirements "mem", "jlabel" job label
##################
bsub_cmd = function(cmd, queue, jname = NULL, jlabel=NULL, jgroup = NULL, mem=NULL, group = "cgafolk", cwd = NULL, mc.cores = NULL, deadline = F)
  {
    if (is.null(jname) & is.null(names(cmd)))
      jname = 'job'

    if (length(jname) != length(cmd))
      jname = rep(jname, length(cmd))
    
    if (!is.null(jname))
      names(cmd) = dedup(jname)    
                            
    qjname = paste( "\"", names(cmd), "\"", sep="" )
    qjout = paste( "\"", names(cmd), ".bsub.out", "\" ", sep="" )
    qjerr = paste( "\"", names(cmd), ".bsub.err", "\" ", sep="" )
    qjrout = paste( "\"", names(cmd), ".R.out", "\" ", sep="" )
    out_cmd = paste( "bsub -q ", queue, " -o ", qjout, " -e ",  qjerr, " -P ", group);
    if (!is.null(mem)) out_cmd = paste(out_cmd, " -R \"rusage[mem=", mem, "]\" ", sep = "");
    if (!is.null(jlabel)) out_cmd = paste(out_cmd, " -J ", jlabel )
    if (!is.null(jgroup)) out_cmd = paste(out_cmd, " -g ", sub('^\\/*', '/', jgroup))
    if (!is.null(cwd)) out_cmd = paste(out_cmd, " -cwd ", cwd )
    if (!is.null(mc.cores)) out_cmd = paste(out_cmd, sprintf(" -n %d,%d -R 'span[hosts=1]'", mc.cores, mc.cores))
    if (deadline) out_cmd = paste(out_cmd, '-sla DEADLINEsla')
    out_cmd = paste(out_cmd," \"",  cmd, "\"", sep = "")
    names(out_cmd)= names(cmd)
    return(out_cmd)
  }


##############
# query_lsf_out
#
# parses "out" and "err" files of jobs with jname root to identify exit status and error codes of jobs
#
##############
lsf_out_query = query_lsf_out = function(dir = NULL, jname = NULL, detailed = F)
{      
  if (!is.null(dir))
    dir = paste(dir, '/', sep = '')
  else
    dir = ''  
  
  input.jname = jname
  
  jname = gsub('\\.bsub.out', '' , dir(dir, '.bsub.out$'))
  tmp.run = paste(normalizePath(dir), '/', jname, sep = '')
  
  if (length(jname)==0)    
    outs = data.frame(jname = NA,
      out.file = NA,
      err.file = NA,
      run.file = NA,
      exit_flag = NA, term_flag = NA, started = NA, reported = NA, hours_elapsed = NA, max_mem = NA, cpu_time = NA, stringsAsFactors = F)
  else
    {
      outs = data.frame(jname = gsub('\\.R$', '', jname),
        out.file = paste(normalizePath(dir),'/', jname, '.bsub.out', sep = ''),
        err.file = paste(normalizePath(dir), '/', jname, '.bsub.err', sep = ''),
        run.file = ifelse(file.exists(tmp.run), tmp.run, NA),
        exit_flag = NA, term_flag = NA, started = NA, reported = NA, hours_elapsed = NA, max_mem = NA, cpu_time = NA, stringsAsFactors = F);


      fn = paste(dir, jname, '.bsub.out', sep = '')
      fn.ex = file.exists(fn);
      if (!any(fn.ex))
        break
      
      tmp = matrix(unlist(lapply(fn[fn.ex],
        function(x)
        {
          y = readLines(x);
          y = split(y, cumsum(grepl('^Sender', y)))
          y = y[[length(y)]]  ## picks "last" dump from lsf to this out file
          return(c(c(grep('^Exited with', y, value = T), grep('^Successfully completed', y, value = T), '')[1],
                   c(grep('^TERM', y, value = T), '')[1],
                   c(gsub('Started at ', '', grep('^Started at', y, value = T)), '')[1],
                   c(gsub('Results reported at ', '', grep('^Results reported at', y, value = T)), '')[1],
                   c(gsub('[ ]+CPU time[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+CPU time', y, value = T)), '')[1],
                   c(gsub('[ ]+Max Memory[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+Max Memory', y, value = T)), '')[1],
                   c(gsub('[ ]+Max Swap[ ]+\\:[ ]+(.*)[ ]+\\S+', '\\1', grep('^[ ]+Max Swap', y, value = T)), '')[1],
                   c(gsub('[ ]+Max Processes[ ]+\\:[ ]+(.*)\\S*', '\\1', grep('^[ ]+Max Processes', y, value = T)), '')[1],
                   c(gsub('[ ]+Max Threads[ ]+\\:[ ]+(.*)\\S*', '\\1', grep('^[ ]+Max Threads', y, value = T)), '')[1]
                   ))

        })), ncol = 9, byrow = T)
      colnames(tmp) = c('exit.flag', 'term.flag', 'started', 'reported', 'cpu.time', 'max.memory', 'max.swap', 'max.cpu', 'max.thr')
      
      TIME.FORMAT = '%a %b %d %H:%M:%S %Y';   
      outs$exit_flag[fn.ex] = tmp[, 'exit.flag']
      outs$term_flag[fn.ex] = tmp[, 'term.flag']
      outs$started[fn.ex] = as.character(as.POSIXct(strptime(tmp[, 'started'], TIME.FORMAT)))
      outs$reported[fn.ex] = as.character(as.POSIXct(strptime(tmp[, 'reported'], TIME.FORMAT)))  
      outs$hours_elapsed = round(as.numeric((as.POSIXct(outs$reported)-as.POSIXct(outs$started))/60), 2)
      outs$cpu_time[fn.ex] = as.numeric(tmp[, 'cpu.time'])
      outs$max_mem[fn.ex] = as.numeric(tmp[, 'max.memory'])

      if (detailed)
        {
          outs$max_swap[fn.ex] = tmp[, 'max.swap']
          outs$max_processes[fn.ex] = tmp[, 'max.processes']      
          outs$max_threads[fn.ex] = tmp[, 'max.threads']
        }
      outs$success = grepl('Success', outs$exit_flag)
      rownames(outs) = dedup(outs$jname)
    }
  
  if (!is.null(input.jname))
    {
      outs = outs[input.jname, ]
      outs$jname = input.jname
      rownames(outs) = dedup(input.jname)
    }
  
  return(outs)
}


###################
# chunk
#
# takes same input as seq (from, to, by, length.out) and outputs a 2 column matrix of indices 
# corresponding to "chunks"
#
###################
chunk = function(from, to = NULL, by = 1, length.out = NULL)
  {
    if (is.null(to))
      {
        to = from;
        from = 1;
      }

    if (is.null(length.out))
      tmp = c(seq(from = from, to = to, by = by), to + 1)
    else
      tmp = c(seq(from = from, to = to, length.out = length.out), to + 1)
    
    out = floor(cbind(tmp[-length(tmp)], tmp[-1]-1))
    
    return(out)
  }

##################
# chunk_and_dump
#
# Takes data frame and chunks into either "k" chunks or so that each chunk has <= "nrow" rows
# then dumps each chunk to file, returning a vector of resulting file paths
#
# If nrow is specified then "k" is ignored. 
##################
chunk_and_dump = function(df, file.prefix = 'chunk', k = 10, nrow = NULL, sep = "\t", quote = F, row.names = F,
  no.dump = F, # no.dump if you want just filenames
  ...) ## ... passed to write.table
  {
    if (!no.dump)
      system(paste('mkdir -p', gsub('\\/([^\\/]+)$', '\\/', file.prefix)))

    if (!is.null(nrow))
      ix = seq(1, nrow(df), nrow)
    else
      ix = round(seq(1, nrow(df), nrow(df)/k))

    if (ix[length(ix)] != nrow(df))
      ix = c(ix, nrow(df)+1)

    fn = paste(file.prefix, ix[1:(length(ix)-1)], sep = ".")

    if (!no.dump)
      {
        if (length(ix)==1)
          write.table(df[1, ], file = fn, sep = sep, quote = quote, row.names = F, ...)
        else
          sapply(1:(length(ix)-1), function(x) write.table(df[ix[x]:(ix[x+1]-1), ], file = fn[x], sep = sep, quote = quote, row.names = F, ...))
      }
    
    return(fn)
  }
      


##################################
# write.bed
#
# dumps segs df with fields $chr, $start, $end, $label to a bed file with filename fn
# also handles segs with chromStart, chromEnd, startbp / endbp, pos1 / pos2 nomenclature
################################## 
write.bed = function(segs, fn, chr = FALSE, header = T, uncollapse = F)
  {

    if (header)
      writeLines(paste("track name = ", fn, "color = 0,0,255"), con = fn)
    else
      system(paste('rm -f', fn))

    segs = standardize_segs(segs);

    if (is.null(segs$label))
      segs$label= paste('seg', 1:nrow(segs), sep = "");

    if (chr)
      segs$chr = paste('chr', gsub('chr', '', segs$chr), sep = "");

    if (uncollapse)
      segs$pos2[segs$pos1==segs$pos2] = segs$pos2[segs$pos1==segs$pos2] + 1;
    
    write.table(segs[, c('chr', 'pos1', 'pos2', 'label')], file = fn, append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE);        
  }

read.bed = function(bed_fn)
{
  bed = read.delim(bed_fn, skip = 1, strings = F, header = FALSE);  
  colnames(bed)[1:4] = c('chr', 'start', 'end', 'name');  
  bed$chr = gsub('chr', '', bed$chr);
  return(bed)
}

read.gff = function(gff_fn)
{
  gff = read.delim(gff_fn, skip = 2, strings = F, header = FALSE);  
  colnames(gff) = c('chr', 'name', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'group');  
  gff$chr = gsub('chr', '', gff$chr);
  return(gff)
}



##################
## Produces (simple) R code calling function named "func" with args in list "argv", prepending with
## source() call to directories in the vector "sources" if specified.
##
## NOTE: args in ... can be lists or vectors consisting of numerical values or characters.  Lists can have named fields.  
## These will be assigned in a "hard coded" way in the Rcode, so these should be ideally scalars or
## pretty short vectors / lists. 
##
## For code to run properly, the names of "argv" must correspond to argument names of "func", or
## if the list has unnamed fields then they must be ordered in the order of the function args.  
##
## Useful for dumping tmp code files for farming when there are many arguments being passed
##################
func_code = function(func, sources = c(), ...)
  { 
    out = "";

    argv = list(...);
    
    if (length(sources)>0)
      {
        out = sprintf('%s%s\n', out, paste("source(\"", sources, "\")", sep = "", collapse = "\n"));
      }

    argv_strings = vector(mode="character", length = length(argv));

    for (i in 1:length(argv))
      {
        this_arg = eval(argv[[i]]); # need to eval if data frame slice passed down as vector (i.e. as "call")
        
        if (is.list(this_arg) & is.null(dim(this_arg))) # checks we have a bona fide list ie not a data frame
          {
            if (max(unlist(lapply(this_arg, length)))>1)
            {
              print("Error: nested list arguments not allowed in argv");
              return( NA );
            }

            list_strings = as.vector(this_arg)
            
            chars = unlist(lapply(this_arg, is.character));  # put quotes around char list items
            list_strings[chars] = paste('\"', list_strings[chars], '\"', sep = "");
            
            if (!is.null(names(this_arg))) # take care of named list items if exists
              {
                named = names(this_arg) != "";        
                list_strings[named] = paste(names(this_arg)[named], " = ", list_strings[named],  sep = "");  # prepend "name=" to named items of list
              }
            
            argv_strings[[i]] = sprintf("list(%s)", paste(list_strings, collapse = ", "));  # pre-pend list constructor and comma concat list items            
          }
        else if (is.vector(this_arg) & is.null(dim(this_arg))) # make sure we have vector and not an array
          {
            vec_strings = this_arg;
            if (is.character(this_arg))
              vec_strings = paste('\"', vec_strings, '\"', sep = "");
            
            if (length(vec_strings)>1) # use c() if we have a vector greater than length 1
              argv_strings[i] = sprintf("c(%s)", paste(vec_strings, collapse = ","))
            else
              argv_strings[i] = vec_strings;
          }
        else if (is.null(this_arg))
          argv_strings[i] = 'NULL'            
        else          
          {
            print("Error: unsupported data type in argv");
            return( NA );
          }
      }
    
    if (!is.null(names(argv))) # take care of named args if exist
      {
        named = names(argv) != "";  
        argv_strings[named] = paste(names(argv)[named], " = ", argv_strings[named], sep = ""); 
      }

    out = sprintf('%s\n%s(%s)\n', out, func, paste(argv_strings, collapse = ",\n "));
    out
  }

#######################
# first.bigger
# 
# Given two <ordered> vectors a and b, returns vector v such that v[i] = index of first element in b that is greater than a[i]
#
#######################
first.bigger = function(a, b)
  {
    j = 0;
    out = vector(mode = "numeric", length = length(a))+NA;

    for (i in 1:length(a))
      {
        while (j<length(b) & b[j+1] <= a[i]) j=j+1;
        if (j < length(b)) out[i] = j+1;
      }
    return(out)
  }

uind = function(list)
  {
    uel = unique(list)
    sapply(uel, function(x){list(which(list==x))})
  }

rand = function(x, y=1)
{
  if (is.numeric(x))
    {
      vdim = x*y;
      array(runif(vdim), dim = c(x, y))	
    }
  else
    {
      ## is a matrix
      array(runif(prod(dim(x))), dim = dim(x))
    }
}

#############################
# readAmpDel = function(fn)
#
# reads GISTIC peak output
############################
readAmpDel = function(fn){
	lines <- readLines(fn, n = 10)
	items = noquote(strsplit(lines[grep("wide peak boundaries", lines)], '\t')[[1]])
	parsed = lapply(items[2:length(items)], function(x){ gsub('^chr(..?):([0-9]+)-([0-9]+)$', '\\1\t\\2\t\\3', x, perl = TRUE)})
	parsed2 = rapply(parsed[lapply(parsed, length)>0] , function(x){strsplit(x, '\t')}[[1]]);
	CNAS = data.frame(t(matrix(parsed2, ncol = length(parsed2)/3, nrow = 3)), stringsAsFactors = FALSE)
	names(CNAS) = c('chr', 'startbp', 'endbp')
	CNAS$startbp = as.numeric(CNAS$startbp)
	CNAS$endbp = as.numeric(CNAS$endbp)
	items = strsplit(lines[grep("cytoband", lines)], '\t')[[1]]
        CNAS$name = items[2:length(items)];
        items = strsplit(lines[grep("q value", lines)], '\t')[[1]]
        CNAS$q.val = as.numeric(items[2:length(items)]);
        items = strsplit(lines[grep("residual q value", lines)], '\t')[[1]]
        CNAS$res.q.val = as.numeric(items[2:length(itemips)]);        
	CNAS
}

#############################
# rrbind = function(df1, df2, [df3 ... etc], )
#
# like rbind, but takes the intersecting columns of the dfs
#
# if union flag is used then will take union of columns (and put NA's for columns of df1 not in df2 and vice versa)
############################
rrbind = function(..., union = T)
  {
    dfs = list(...);  # gets list of data frames
    dfs = dfs[!sapply(dfs, is.null)]
    dfs = dfs[sapply(dfs, ncol)>0]
    names.list = lapply(dfs, names);
    cols = unique(unlist(names.list));
    unshared = lapply(names.list, function(x) setdiff(cols, x));
    expanded.dfs = lapply(1:length(dfs), function(x) {dfs[[x]][, unshared[[x]]] = NA; return(dfs[[x]])})
    out = do.call('rbind', expanded.dfs);

    if (!union)
      {
        shared = setdiff(cols, unique(unlist(unshared)))
        out = out[, shared];
      }
     
   return(out)
  }

#############################
# read.delim.cat(paths, ... [read.delim arguments])
#
# takes a vector of tab delimited file paths and concatenates them into a
# single data frame (takin union of identically named / numbered columns as a default)
############################
read.delim.cat = function(paths, skip = NULL, cols = NULL, include.paths = T, cores = NULL, ...)	
  {
    if (is.null(skip))
      skip = rep(0, length(paths))
    
    paths[is.na(paths)] = ""; 
    does.not.exist = !file.exists(paths);
    
    if (any(does.not.exist))
      warning(sprintf('Ignoring %s paths that do not exist on the file system.', length(which(does.not.exist))))
    paths = paths[!does.not.exist];
    
    if (length(skip) ==1)
      skip = rep(skip, length(paths))
    else
      skip = skip[!does.not.exist]
    
    # scope out files to filter out those with 0 rows and find common columns
    if (!is.null(cores))
      dfs = mclapply(1:length(paths),
        function(x) {tmp.df = read.delim(paths[x], skip = skip[x], ...);
                     if (nrow(tmp.df) != 0) cbind(data.frame(source.path = paths[x]), tmp.df) else data.frame() }, mc.cores = cores)
    else
      dfs = lapply(1:length(paths),
        function(x) {tmp.df = read.delim(paths[x], skip = skip[x], ...);
                     if (nrow(tmp.df) != 0) cbind(data.frame(source.path = paths[x]), tmp.df) else data.frame() })


    dfs = dfs[sapply(dfs, nrow)!=0];    
    
    if (length(dfs)==0)
      return(NULL);
   
    out = do.call('rrbind', dfs)

    if (!is.null(cols))
      out = cbind(out[,1], out[, cols]);

    if (include.paths)
      names(out)[1] = 'source.path'
    else
      out = out[,-1]
   
    return(out)
  }


####################
# Retrieves all snps in r2 "ldthresh" of snps in region defined by data.frame region with fields chr pos1 (and optional pos2) fields
# 
# Example:
#
# ldQuery(data.frame(chr ='chr1', pos1 = 32323232, pos2 = 32343232), 0.9);
#
#
# RMysql installation cmd
#  install.packages('RMySQL', type='source', lib = '~/R/x86_64-unknown-linux-gnu-library/2.8/', configure.args = "--with-mysql-inc='/local/util/include/mysql/' --with-mysql-lib='/local/util/lib/mysql/'")
####################
ldQuery = function(region = NULL, ldthresh)
{
  require('RMySQL')
  con = dbConnect(MySQL(),
    user = "hapmap", password = "qwerasdf", dbname="hapmap", host=MARCIN.COMPUTER)

  if (is.null(region$pos2))
      region$pos2 = region$pos1;

  if (!any(grep('chr', region$chr)))
    region$chr = paste('chr', region$chr, sep = "");

  q1 = sprintf('select * from ld_CEU where chr = "%s" and r2 > %f and pos1 between %s and %s', 
    region$chr, ldthresh, region$pos1, region$pos2)
  
  q2 = sprintf('select * from ld_CEU where chr = "%s" and r2 > %f and pos2 between %s and %s;', 
    region$chr, ldthresh, region$pos1, region$pos2)
  
  out = rbind(dbGetQuery(con, q1), dbGetQuery(con, q2))

  dbDisconnect(con)

  return(out[!duplicated(out),])  
}	

####################
# Retrieves all pairwise CEU ld values for snps in list snp
#
# ie ldSNP(snps)
#
####################
ldSNP = function(snps)
{
  require('RMySQL')
  con = dbConnect(MySQL(),
    user = "hapmap", password = "qwerasdf", dbname="hapmap", host="vmfe7-85b")

  q = sprintf('select * from ld_CEU where snp1 in (\"%s\") and snp2 in (\"%s\")', 
     paste(snps, collapse = '","'), paste(snps, collapse = '","'))

print(q)
  
  out = dbGetQuery(con, q);

  dbDisconnect(con)

  return(out)
}

####################
# ldSNP = function(snps, pop)
## Retrieves mafs for snp vector "snps" in HapMap pop "pop"
#
####################
af_SNP = function(snps, pop = 'CEU')
{
  require('RMySQL')
  con = dbConnect(MySQL(),
    user = "hapmap", password = "qwerasdf", dbname="hapmap", host="vmfe7-85b")

  q = sprintf('select * from af where snp in (\"%s\") and pop = \"%s\"', 
     paste(snps, collapse = '","'), pop)
  
  out = dbGetQuery(con, q);
  rownames(out) = out$snp;
  out = out[snps,];
  out$snp = snps;
  rownames(out) = NULL;
  
  dbDisconnect(con)

  return(out)
}

###############################
# fisher.plot
#
# Plots fisher contingency table 
#
###############################
fisher.plot = function(O)
  {
    require(plotrix)
    fish = fisher.test(O)
    plot.new();
    par(usr = c(0, 1, 0, 1));
    addtable2plot(0,0.5, O, display.colnames = T, display.rownames = T);
    text(0.5, 0.44, sprintf(paste('P = %0.', floor(-log10(fish$p.value))+2, 'f\nOR = %0.2f [%0.2f-%0.2f]', sep  = ""), fish$p.value, fish$estimate, fish$conf.int[[1]], fish$conf.int[[2]]))
  }

###############################
# fisher.combined
#
# Computes fisher combined p value for a matrix of p values where the columns correspond to individual tests
# rows correspond to hypotheses.
#
###############################
fisher.combined = function(Ps)
{
  if (is.vector(Ps))
    return(Ps)

  return(pchisq(rowSums(-2*log(Ps)), 2*ncol(Ps), lower.tail = F))
}

####################
# retrieves the region from CEU hapmap with r2>ld to the query region
####################
ldHood = function(region, ld)
  {
    if (is.null(region$pos2))
       {
         region$pos2 = region$pos1;
       }

     temp = ldQuery(region, ld);
     
     if (dim(temp)[1]>0)
       {
         pos = rbind(temp$pos1, temp$pos2);         
         out = data.frame(chr=temp$chr[1], pos1=min(pos), pos2=max(pos))
       }
     else
       {
         out = data.frame(chr=region$chr, pos1=region$pos1, pos2=region$pos2)
       }

    return(out)
  }

lsf.status = function()
  {
    TEMP.FILE = 'tmp.432124kiijoij1oij2oii1ijoi202082892h4iojh12oij12';
#    system(paste("bjobs -u all | awk \'{print", paste('$', c(1:7), sep = "", collapse = "\"\t\""), "\"\t\"", paste('$', c(8:10), sep = "", collapse = "\" \""), "}\' > ", TEMP.FILE, sep = ""))
    system(paste("bjobs -u all > ", TEMP.FILE, sep = ""))
        
    header = readLines(TEMP.FILE,1);
    headers = strsplit(header[[1]], "\\s+");
    headers = sub("\\_", "\\.", headers[[1]]);
    
    w = get.field.widths(header[[1]]);
    w[length(w)] = 1000;
   
#    bout = read.delim(TEMP.FILE, stringsAsFactors = FALSE, header = FALSE, skip = 1, col.names = c('JOBID','USER','STAT','QUEUE','FROM_HOST','EXEC_HOST','JOB_NAME','SUBMIT_TIME'));
    bout = read.fwf(TEMP.FILE, w, stringsAsFactors = FALSE)
    system(paste('rm ', TEMP.FILE));

    for (j in 1:dim(bout)[2])
        bout[,j] = trim(bout[,j]);       

    bout$SUBMIT.TIME = as.POSIXct(paste(bout$SUBMIT.TIME,  format(Sys.Date(), "%Y")), format =  '%b %d %H:%M %Y');
    
    bout
  }

lsf.usersum = function(lsfstat, attr = "STAT")
  {
    tab = table(lsfstat$USER, lsfstat[, attr]);
    tab = tab[order(-rowSums(tab)),]
    tab
  }

get.fwf.widths = function(file, skip=0)
  {
    l = readLines(file,skip+1);

    w = get.field.widths(l[[length(l)]]);    
  }


############################
# dev.all.off
#
# kills all windows
#
#############################
dev.all.off = function()
  {
    sapply(dev.list(), dev.off)
  }


# given logical vector returns data frame of starts and end indices of TRUE run lengths
run.lengths = function(vec)
  {
    out = data.frame();
    
    vec = c(vec, FALSE)
    last.bit = FALSE;
    this.bit = FALSE;
    last.true = as.double(this.bit);
    
    for (i in 1:length(vec))
      {
        last.bit = this.bit;
        this.bit = vec[i];
        
        if (!this.bit & last.bit)
          {
            out = rbind(out, data.frame(start = last.true, end = i-1))
          }
        
        if (this.bit & !last.bit)
          {
            last.true = i;
          }
      }
    
    out
  }

get.field.widths = function(str)
  {
    spl = strsplit(str, "")
    non.space = is.na(match(spl[[1]], " "));

    runs = run.lengths(non.space);

    out = c();
    
    if (dim(runs)[1]==1)
      {
        out = length(spl[[1]]);
      }
    else if (dim(runs)[1]>1)
      {         
        out = c(runs$start[2:dim(runs)[1]], length(spl[[1]])+1) -
          c(1, runs$start[2:dim(runs)[1]]);
      }

    out
  }

trim = function(str)
  {
    nms = names(str)
    str = gsub("^\\s+", "", str);
    str = gsub("\\s+$", "", str);
    return(structure(str, names = nms))
  }


chiblocks = function(b1, b2)
{
	tally = matrix(0, nrow=dim(hapblocks)[1], ncol =2); 
	tally[match(b1, hapblocks$BLOCKNAME), 1] = 1; 
	tally[match(b2, hapblocks$BLOCKNAME), 2] = 1;
	tab = table(as.factor(tally[,1]), as.factor(tally[,2])); 
	chisq.test(tab)
}

linblocks = function(b1, b2)
{
	tally = data.frame(matrix(0, nrow=dim(hapblocks)[1], ncol =2))
	names(tally) = c('b1', 'b2')
	
	ind = match(b1, hapblocks$BLOCKNAME)
	tally[ind[!is.na(ind)], 'b1'] = 1 

	ind = match(b2, hapblocks$BLOCKNAME)
	tally[ind[!is.na(ind)], 'b2'] = 1 

	lm('b2 ~ b1', tally)
}


listind = function(query, thelist)
{
     out = vector("list", length(query))
     for (i in (1:length(thelist)))
       {
          m = which(!is.na(match(query, thelist[[i]])))
          if (length(m)>0) {
            for (j in (1:length(m)))
              {
                out[[m[j]]] = c(out[[m[j]]], i)
              }
          }
       }
     out
}

#################
# write.tab - writes tab delimited no quotes without row names table (passes remaining arguments to write.table)
#
# = write.table(sep = "\t", quote = F, row.names = F)
#################
write.tab = function(..., sep = "\t", quote = F, row.names = F)
  {
    write.table(..., sep = sep, quote = quote, row.names = row.names)
  }


#################
# write.htab - writes data frame to pretty formatted rtable
#
#
#################
write.htab = function(tab, file = NULL,
  title = NULL, # text to be written in bold above the table  
  footer = NULL, # text to be writen in bold below the table
  highlight = NULL,  #vector of row indices of the table to highlight
  row.names = TRUE,  # includes row labels
  col.names = TRUE, # includes col labels
  high.color = 'yellow', # highlight color to use 
  row.colors = c('lightgray', 'white'), # alternating colors to shade data rows  
  header.colors = c('#4A4A4A', 'white'), # two element vector specifying background and text colors for header row, respectively,
  data.size = 15, # font size in px for data, title, and footer
  title.size = 15, footer.size = 20, header.size = round(1.1*data.size))
  {    
    require(hwriter)
    require(gplots)
    
    if (class(tab) != 'data.frame')
      tab = as.data.frame(tab)

    if (is.null(rownames(tab)))
      row.names = F;

    if (!is.null(file))
      {
        if (!grepl('($~)|($\\/)', file))
          file = paste('~/public_html/', file, sep = '')
      }
    else
      file = DEFAULT.HTAB.FILE    
    
    tab[is.na(tab)] = '';
    tab = tab[1:nrow(tab), , drop = FALSE];  #not sure why this is necessary, but deflects occasional weird R bug

    if (any(lix <<- sapply(names(tab), function(x) is.list(tab[, x]))))
      for (i in which(lix))
        tab[, i] = sapply(tab[, i], function(x) paste(x, collapse = ','))
    
    dir.create(dirname(normalizePath(file.dir(file))), recursive=T, showWarnings = F)
    p = openPage(file, link.css = 'hwriter.css')
    if (!is.null(title))
      hwrite(title, p, style = sprintf('font-weight:bold; font-size:%spx; margin-top;50px', title.size), center = TRUE, div = TRUE, br = TRUE);

    row.bgcolor = as.list(as.character(col2hex(row.colors)[(1:nrow(tab))%%length(row.colors)+1]));
    names(row.bgcolor) = rownames(tab)
    if (!is.null(highlight))
      row.bgcolor[rownames(tab[highlight,, drop = FALSE])] = list(col2hex(high.color));

    row.bgcolor = c(col2hex(header.colors[1]), row.bgcolor)

#    if (row.names)
      col.bgcolor = col2hex(header.colors[1])
    
    col.style = sprintf('font-weight:bold; font-size:%spx; color:%s; text-align:center', header.size, col2hex(header.colors[2]));
    
    row.style = rep(sprintf('font-size:%spx; text-align:center', data.size), nrow(tab))
    names(row.style) = rownames(tab)
    row.style = c(list(sprintf('font-weight:bold; font-size:%spx; color:%s; text-align:center', header.size, col2hex(header.colors[2]))), row.style)
    
    hwrite(tab, p, row.style = row.style, col.style = col.style, col.bgcolor = col.bgcolor, row.names = row.names, col.names = col.names,
           row.bgcolor = row.bgcolor, table.frame = 'void', table.style = 'margin-left: 30px; margin-top: 30px', br = TRUE)
    if (!is.null(footer))
      hwrite(footer, p, style = sprintf('font-weight:bold; text-align:center; font-size:%spx; margin-top;50px', footer.size), center = TRUE, div = TRUE);
    closePage(p)
  }

###############
# writecols
#
# Takes character vector and dumps into k delimited columns
###############
writeCols = function(v, k = 3, sep = "\t", file = "")
  {
    rows = ceiling(length(v)/k)
    out = matrix(ncol = k, nrow = rows, byrow = TRUE)
    out[1:length(v)] = v;
    out[is.na(out)] = "";
    write.table(as.data.frame(out), file = file, sep = sep, quote = F, row.names = F, col.names = F)
  }


##########################
# col.scale
#
# Assigns rgb colors to numeric data values in vector "x".. maps scalar values
# in val.range (default c(0,1)) to a linear color scale of between col.min (default white)
# and col.max (default black), each which are length 3 vectors or characters.  RGB values are scaled between 0 and 1. 
#
# Values below and above val.min and val.max are mapped to col.max and col.max respectively
##########################
col.scale = function(x, val.range = c(0, 1), col.min = 'white', col.max = 'black', na.col = 'white',
  invert = F # if T flips rgb.min and rgb.max
  )
  {
    if (!is.numeric(col.min))
      if (is.character(col.min))
        col.min = col2rgb(col.min)/255
      else
        error('Color should be either length 3 vector or character')

    if (!is.numeric(col.max))
      if (is.character(col.max))
        col.max = col2rgb(col.max)/255
      else
        error('Color should be either length 3 vector or character')

    col.min = as.numeric(col.min);
    col.max = as.numeric(col.max);

    x = (pmax(val.range[1], pmin(val.range[2], x))-val.range[1])/diff(val.range);
    col.min = pmax(0, pmin(1, col.min))
    col.max = pmax(0, pmin(1, col.max))

    if (invert)
      {
        tmp = col.max
        col.max = col.min
        col.min = tmp
      }

    nna = !is.na(x);
    
    out = rep(na.col, length(x))
    out[nna] = rgb((col.max[1]-col.min[1])*x[nna] + col.min[1],
        (col.max[2]-col.min[2])*x[nna] + col.min[2],
        (col.max[3]-col.min[3])*x[nna] + col.min[3])
    
    return(out)           
  }

###########################
# axis.subplot
#
# function to draw simple horizontal and vertical axes at arbitrary positions on plot with
# decoupled data and plot coordinates (useful when there are many subplots within a single layout item)
#
###########################
axis.subplot = function(side = 1, # can be 1 (horizontal) or 2 (vertical)
  data.range, # (data coordinates) min and max data coordinate value for ends of axis line
              # corresponding to ypos (or xpos) plot coordinates in case of vertical (or horizontal) plot, respectively
  n.ticks = 5, # optimal number of ticks using "pretty" function
  at = NULL, # (data coordinates) tick locations same as in original axis function, over-rides n.ticks 
  ypos = NULL, # (plot coordinates) single value for horizontal, two values (range) for vertical correponding to ypos of axis line
  xpos = NULL, # (plot coordinates) single value for vertical, two values (range) for horizontal corresponding to xpos of axis line
  lwd = 1,  
  lwd.ticks = lwd,
  col = 'black',
  col.ticks = col,
  tick.length = 0.1, # in plot coordinates
  tick.orient = 1, # if >0 then facing "positive" direction, otherwise "negative"
  srt.tick = 0, # rotation of tick.labels relative to "default" (which is 0 for horizontal and 90 for vertical axis)
  cex.tick = 1, # size of tick labels
  adj.tick = c(0.5, 0.5), # justification of tick.labels (overrides defaults)
  ...) # arguments to line)
{  
  if (side == 1 & !(length(xpos)==2 & length(ypos)==1))
    stop('ypos must be of length 1 and xpos of length 2, for side = 1')
  else if (side == 2 & !(length(xpos)==1 & length(ypos)==2))
    stop('ypos must be of length 2 and xpos of length 1, for side = 2')
  else if (!(side %in% c(1, 2)))
    stop('side must be = 1 or = 2')
  
  axis.lim = data.frame(x = xpos, y = ypos);

  if (is.null(at))
    at = pretty(data.range, n = n.ticks)

  # draw axis line
  lines(axis.lim$x, axis.lim$y, lwd = lwd, col = col, ...)

  # draw ticks
  if (tick.length <= 0)
    tick.length = -tick.length;
  
  if (side == 1)
    {
      at.plot = coord.subplot(at, data.range = data.range, plot.range = xpos)
      segments(at.plot, rep(ypos, length(at)), y1 = rep(ypos, length(at))+tick.length, lwd = lwd, col = col, ...)
      text(rep(at.plot, 2), rep(ypos, length(at))+1.5*tick.length, at, adj = adj.tick, cex = cex.tick, srt = 0+srt.tick)
    }
  else if (side == 2)
    {      
      at.plot = coord.subplot(at, data.range = data.range, plot.range = ypos)
      segments(rep(xpos, length(at)), at.plot, x1 = rep(xpos, length(at))+tick.length, lwd = lwd, col = col, ...)
      text(rep(xpos, length(at))+1.5*tick.length, rep(at.plot, 2), at, adj = adj.tick, cex = cex.tick, srt = -90+srt.tick)      
    }      
}

##############################
# coord.subplot
#
# Scales data coordinates in given "data.range" to plot coordinates in "plot.range".  Can be used for x or y coordinates.
# Data coordinates that map outside of the data.range are given NA values. 
#
# Useful for doing many subplots on a single plot axis. 
# e.g. if we want to plot genomic positions 54,100,000 to 56,200,000 on plot coordinates 0.1 to 0.4
##############################
coord.subplot = function(data, # any length data vector  
              data.range, # length 2 vector 
              plot.range # length 2 vector
              )
{
  data[data<data.range[1] | data>data.range[2]] = NA;
  return(((data-data.range[1])/diff(data.range))*diff(plot.range)+plot.range[1])
}

############################
# strwrap2
#
# Wrapper around strwrap that returns a vector of characters the same length as input with new lines
# inserted in word breaks to obey a provided string width.
#
#############################
strwrap2 = function(x, width, sep = NULL, newline = '\n', ...)
  {
    return(sapply(strwrap(x, width = width, simplify = F),
                  function(x) paste(x, collapse = newline)))
  }

############################
# Capitalize first letter of each character element of vector "string"
# 
##############################
capitalize = function(string, un = FALSE)
{
    if (!un)
      {
        capped <- grep("^[^A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- toupper(substr(string[capped],1, 1))
      }
    else
      {
        capped <- grep("^[A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- tolower(substr(string[capped],1, 1))
      }
    
    return(string)
}

############################
# Performs fisher test on cols of matrix / df x vs cols of matrix / df y
#
# returns list with ncol(x) by ncol(y) matrices $p and $or denoting the p value and odds ratio of the result of the
# fisher test on col i of x and col j of y
#
# If y is not provided, will correlate rows of x with themselves. 
#############################
fisher.pairwise = function(x, y = x)
  {
    p = or = matrix(NA, nrow = ncol(x), ncol = ncol(y), dimnames = list(colnames(x), colnames(y)))

    if (nrow(x) != nrow(y))
      stop('x and y must have the same number of rows')
    
    logical.x = which(sapply(1:ncol(x), function(i) is.logical(x[,i])))
    logical.y = which(sapply(1:ncol(y), function(i) is.logical(y[,i])))

    if (length(logical.x)==0 | length(logical.y)==0)
      warning('No logical columns found')

    for (i in logical.x)
          for (j in logical.y)
            {
              O = table(x[,i], y[,j])
              if (min(dim(O))>1)
                {
                  res = fisher.test(O)          
                  or[i, j] = res$estimate
                  p[i, j] = res$p.value
                }
            }

    out = list(p = p, or = or)
    return(out)
  }

#################
# lighten
#
# lightens / darkens colors by brighness factor f (in -255 .. 255) that will make lighter if > 0 and darker < 0
################
lighten = function(col, f)
  {
    M = col2rgb(col)
    return(apply(matrix(pmax(0, pmin(255, M + f*matrix(rep(1, length(M)), nrow = nrow(M)))), ncol = length(col))/255, 2, function(x) rgb(x[1], x[2], x[3])))
  }


################
# fin.lines
#
# returns par('fin') in terms of lines - helps margin computation
#
################
fin.lines = function() par('fin')*(par('mar')/par('mai'))[1];  # help us make sure we don't encroach on figure margins


################
# strsplit2
#
# Strsplit when there are two layers of separators (sep1, sep2) and one needs to extract
# a collapsed vector of subitem j for all items i.
#
# Takes in a character vector and outputs a list of "separated" items
################
strsplit2 = function(x, sep1 = ",", sep2 = " ", j = 1)
  {
    return(lapply(strsplit(x, sep1), function(y) sapply(strsplit(y, sep2), function(z) z[[j]])))
  }


################ 
# timestamp
#
################
timestamp = function()
  {
    return(gsub('[^\\d]', '', as.character(Sys.time()), perl = T))
  }


###############
# img.link
#
# Returns vector of html image links to files "file" with text "caption"
#
# if embed = T, then will make img link, and additional arguments may be supplied to image tag (eg height, width)
################
img.link = function(file, caption = NULL, embed = F, ...) {
	if (is.null(caption)) {
		caption = ''
	}
	
	if (!embed) {
		return(paste('<a href = \"', file, '\">', caption, '</a>'))
	} else {
		args = list(...);
		if (length(args)>0) {
			to.quote = is.numeric(unlist(args))+1
			more.args = paste((paste(', ', names(args), " = ", c('\"', '')[to.quote], unlist(args), c('\"', '')[to.quote], sep = "")), collapse="")
		} else {
			more.args = ''
		}
		
		parts = strsplit(basename(file), "\\.")[[1]]
		file.ext = parts[length(parts)]
		if (file.ext == "tif") { # handles tif images which won't display in chrome with a regular img src tag
			out = paste('<embed src = \"', file, '\", type = "image/tiff"', more.args, ' negative=yes>', sep = "")
		} else {
			out = paste('<img src = \"', file, '\", alt = \"', caption, '\"', more.args, '>', sep = "")
		}
		
		return(out)
	}
}


#########
# html.link
#
# returns text with html link
#########
html.link = function(href, text = NULL)
  {
    if (is.null(text))
      text = file.name(href)
    
    return(mapply(function(x,y) html.tag('a', href = x,  text = y), href, text))
  }


#########
# html.tag
#
# makes a open and close html tag with optional text inside and optional (named) vector of name=value pairs contained inside of opening tag
#########
html.tag = function(tag, text = NULL, collapse = '\n', ## collapse is what will separate the tag and each "line" of text, alternate can be " "
  ...)
  {
    flags = unlist(list(...))
      
    if (!is.null(flags))
      {
        if (is.null(names(flags)))
          flag.str = ""
        else
          flag.str = paste(" ", paste(paste(names(flags), paste('"', flags, '"', sep = ""), sep = "=")), collapse = "")
      }
    else
      flag.str = ""
    
    return(paste(sprintf('<%s%s>', tag, flag.str), paste(text, collapse = collapse), sprintf('</%s>', tag), sep = collapse))
  }


################
# set.comp
#
# Compares two sets and outputs data frame with "left", "middle", "right" members
################
set.comp = function(s1, s2)
  {
    universe = sort(union(s1, s2));
    out = data.frame(stringsAsFactors = F)       
    tmp.comp = list(left = sort(setdiff(s1, s2)), middle = sort(intersect(s1, s2)), right = sort(setdiff(s2, s1)))
    lapply(names(tmp.comp), function(x) out[1:length(tmp.comp[[x]]),x] <<- tmp.comp[[x]])
    out[is.na(out)] = ''
    return(out)
  }


################################
# dedup
#
# relabels duplicates in a character vector with .1, .2, .3
# (where "." can be replaced by any user specified suffix)
################################
dedup = function(x, suffix = '.')
{
  dup = duplicated(x);
  udup = unique(x[dup])
  udup.ix = lapply(udup, function(y) which(x==y));
  udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
  out = x;
  out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
  return(out)  
}

########################
# install.packages.bioc
#
# shortcut to install bioconductor packages
########################
install.packages.bioc = function(pkg)
  {
    source('http://bioconductor.org/biocLite.R')
    sapply(pkg, biocLite)
  }

##########################
# install.packages.github
#
# shortcut to install github packages
##########################
install.packages.github = function(pkg, username, branch)
  {
    library('devtools')
    install_github(repo = pkg, username = username, branch = branch)    
  }

#############################
# levapply
#
# Applies FUN locally to levels of x and returns vector of length()
# (eg can do a "local" order within levels)
#
#############################
levapply = function(x, by, FUN = 'order')
  {
    if (!is.list(by))
      by = list(by)
    
    f = factor(do.call('paste', c(list(sep = '|'), by)))
    ixl = split(1:length(x), f);
    ixv = lapply(ixl, function(y) x[y])    
    res = structure(unlist(lapply(ixv, FUN)), names = unlist(ixl))
    out = rep(NA, length(x))
    out[as.numeric(names(res))] = res;
    return(out)
  }


####################
# cytoscape
#
# shortcut to connect to cytoscape running on Marcin's mac via RPC
#
# graph must be  igraph, or adjacency matrixx
####################
cytoscape = function(graph = NULL, sessionName = 'M-ski', host = MARCIN.COMPUTER, port = 9000, display = T, layout = 'degree-circle', verbose = T, ...)
  {
    require(RCytoscape)
    require(igraph)
    

    if (is(graph, 'matrix'))
      graph = graph.adjacency(graph, weighted = 'weight');
    
    if (is(graph, 'igraph'))
      {
        if (!is.null(E(graph)$weight))
          E(graph)$weight = 1

        if (!is.null(E(graph)$arrow.shape))
          E(graph)$arrow.shape = 'ARROW'
        
        graph.nel = igraph2graph(graph)
      }
    
    cw = new.CytoscapeWindow(sessionName, host = host, rpcPort = port, graph = graph.nel,  ...)
    
    if (display)
      {
        displayGraph(cw)
        setDefaultBackgroundColor(cw, col2hex('white'))
        
        eG = paste(get.edgelist(graph)[,1], get.edgelist(graph)[,2], sep = '~')
        ceG = cy2.edge.names(cw@graph)

        if (verbose)
          cat('Setting line styles\n')
        
        if ('line.style' %in% list.edge.attributes(graph))
          {
            uls = setdiff(E(graph)$line.style, NA)
            setEdgeLineStyleRule(cw, 'line.style', uls, uls)
          }

        if (verbose)
          cat('Setting arrow shape\n')
        
        if (is.directed(graph))
          if ('arrow.shape' %in% list.edge.attributes(graph))
            setEdgeTargetArrowRule(cw, 'arrow.shape', unique(E(graph)$arrow.shape), unique(E(graph)$arrow.shape))

        if (verbose)
          cat('Setting edge color\n')
        
        if ('col' %in% list.edge.attributes(graph))
          {
            uc = setdiff(unique(E(graph)$col), NA);
            setEdgeColorRule(cw, 'col', uc, uc, mode = 'lookup')
          }

        if (verbose)
          cat('Setting edge width\n')
        
        if ('width' %in% list.edge.attributes(graph))
          {
            uw = setdiff(E(graph)$width, NA)
            setEdgeLineWidthRule(cw, 'width', as.character(uw), uw)
          }

        if (verbose)
          cat('Setting node size\n')
        
        if ('size' %in% list.vertex.attributes(graph))
          {
            us = setdiff(unique(V(graph)$size), NA)
            setNodeSizeRule(cw, 'size', us, us, mode = 'lookup')
          }

        if (verbose)
          cat('Setting node color\n')
        
        if ('col' %in% list.vertex.attributes(graph))
          {
            uc = setdiff(unique(V(graph)$col), NA)
            setNodeColorRule(cw, 'col', uc, uc, mode = 'lookup')
          }

        if (verbose)
          cat('Setting node labels\n')
        
        if ('label' %in% list.vertex.attributes(graph))
          setNodeLabelRule(cw, 'label')

        if (verbose)
          cat('Setting node shapes\n')
        
        if ('shape' %in% list.vertex.attributes(graph))
          {
            us = setdiff(unique(V(graph)$shape), NA)
            setNodeShapeRule(cw, 'shape', us, us, default = 'ELLIPSE')
          }
        
        if (verbose)
          cat('Setting node width\n')
        
        if ('border.width' %in% list.vertex.attributes(graph))
          {
            ubw = setdiff(unique(V(graph)$border.width), NA)
            setNodeBorderWidthRule(cw, 'border.width', ubw, ubw)
          }
        
        if (all(c('x', 'y') %in% list.vertex.attributes(graph)))
          {
            good.ix = !is.na(V(graph)$x) & !is.na(V(graph)$y)
            if (any(good.ix))
              setNodePosition(cw, V(graph)$name[good.ix], V(graph)$x[good.ix], V(graph)$y[good.ix])
          }
        else            
          layoutNetwork(cw, layout)

        
        redraw(cw)        
      }

    return(cw)
  }


####################
# igraph2graph
#
# Converts igraph object into janky graphNEL object (for visualization in cytoscape)
# and populates all edge features both via the edgeL and as NodeAttributes for visualization
# in cytoscape
####################
igraph2graph = function(g)
  {
    require(igraph)
    require(graph)
    require(RCytoscape)
    
    if (class(V) != 'function' | class(E) != 'function')
      stop('Namespace conflict - either V() or E() no longer mapping to igraph functions')

    if (!is.null(V(g)$name))
      node.labels = V(g)$name
    else
      node.labels = as.character(V(g));
    
    edge.df = structure(as.data.frame(get.edgelist(g), stringsAsFactors = F), names = c('vertices', 'edges'))
    if (length(list.edge.attributes(g))>0)
      {
        tmp = do.call('cbind', lapply(list.edge.attributes(g), function(x) get.edge.attribute(g, x)))
        colnames(tmp) = list.edge.attributes(g)        
        edge.df = cbind(edge.df, as.data.frame(tmp, stringsAsFactors = F))        
      }
    if (!is.null(edge.df$weights))
      edge.df$weights = as.numeric(edge.df$weights)
    edge.df[is.na(edge.df)] = "NA"
    edge.df[,1] = as.character(edge.df[,1])
    edge.df[,2] = as.character(edge.df[,2])

    vertex.df = data.frame(vertices = node.labels, stringsAsFactors = F)
    if (length(list.vertex.attributes(g))>0)
      {
        tmp = do.call('cbind', lapply(list.vertex.attributes(g), function(x) get.vertex.attribute(g, x)))
        colnames(tmp) = list.vertex.attributes(g)        
        vertex.df = cbind(vertex.df, as.data.frame(tmp, stringsAsFactors = F))        
      }
        
    ## have to reciprocate edges in undirected otherwise graphNEL will barf
    if (!is.directed(g))
      {
        edge.df.rev = edge.df;
        tmp.col = edge.df[,2]
        edge.df.rev[,2] = edge.df.rev[,1]
        edge.df.rev[,1] = tmp.col;
        edge.df = rbind(edge.df, edge.df.rev)        
      }
    
    edgeL = lapply(split(edge.df, edge.df$vertices), as.list)[node.labels]
    names(edgeL) = node.labels;

    ## retarded GraphNEL object format necessitates these gymnastics
    null.vert = sapply(edgeL, is.null)    
    blank.edge.item = list(structure(rep(list(c()), length(list.edge.attributes(g))+1), names = c('edges', list.edge.attributes(g))))
    edgeL[null.vert] = blank.edge.item
    edgemode = c('undirected', 'directed')[1 + is.directed(g)]
    
    out.g = new('graphNEL', node.labels, edgeL, edgemode = edgemode)
    
    ## populate edge and node attribute for RCytoscape to access
    if (ncol(edge.df)>2)
      {
        attr.cols = names(edge.df)[3:ncol(edge.df)]
        for (attr in attr.cols)          
          {
            if (is.numeric(edge.df[, attr]))                   
              out.g = initEdgeAttribute(out.g, attr, 'numeric', NA)
            else if (is.integer(edge.df[, attr]))
              out.g = initEdgeAttribute(out.g, attr, 'integer', NA)
            else
              {
                cast.numeric = suppressWarnings(as.numeric(edge.df[, attr]))
                cast.character = edge.df[, attr]

                if (any(is.na(cast.numeric) & cast.character != "NA"))
                  {
                    edge.df[, attr] = as.character(edge.df[, attr])
                    out.g = initEdgeAttribute(out.g, attr, 'char', '')
                  }
                else
                  {
                    cast.numeric[is.na(cast.numeric)] = 0
                    edge.df[, attr] = cast.numeric
                    out.g = initEdgeAttribute(out.g, attr, 'numeric', '')
                  }
              }
            edgeData(out.g, edge.df[,1], edge.df[,2], attr) = edge.df[,attr]
          }
      }
        
    if (ncol(vertex.df)>1)
      {
        attr.cols = names(vertex.df)[2:ncol(vertex.df)]
        for (attr in attr.cols)
          {
            if (is.numeric(vertex.df[, attr]))                 
              out.g = initNodeAttribute(out.g, attr, 'numeric', NA)
            else if (is.integer(vertex.df[, attr]))
              out.g = initNodeAttribute(out.g, attr, 'integer', NA)
            else
              {
                cast.numeric = suppressWarnings(as.numeric(vertex.df[, attr]))
                cast.character = suppressWarnings(as.character(vertex.df[, attr]))
                if (any(setdiff(is.na(cast.numeric) & cast.character != "NA", NA)))
                  {
                    vertex.df[, attr] = as.character(vertex.df[, attr])
                    out.g = initNodeAttribute(out.g, attr, 'char', '')
                  }
                else
                  {
                    cast.numeric[is.na(cast.numeric)] = 0
                    vertex.df[, attr] = cast.numeric
                    out.g = initNodeAttribute(out.g, attr, 'char', '')
                  }
              }
            nodeData(out.g, vertex.df[,1], attr) = vertex.df[,attr]
          }
      }
    
    return(out.g)
  }


################################
# cyto2igraph
#
# Converts graph in cytoscape window "cw" into igraph object
################################
cyto2igraph = function(cw)
{
  node.attr = getAllNodeAttributes(cw)
  edge.attr = getAllEdgeAttributes(cw)
  directed = edgemode(cw@graph) == 'directed'

  edge.df = cbind(from = edge.attr$source, to = edge.attr$target,
    edge.attr[, setdiff(names(edge.attr), c('source', 'target'))])

  if ('name' %in% node.attr)
    node.df = node.attr[, c('name', setdiff(names(node.attr), 'name'))]
  else
    node.df = cbind(name = rownames(node.attr), node.attr)

  return(graph.data.frame(edge.df, directed = directed, vertices = node.df))
}

###############################
# brewer.master
#
# Makes a lot of brewer colors using entire brewer palette
#
###############################
brewer.master = function(n, palette = 'Accent')
{
  library(RColorBrewer)
  palettes = list(
    sequential = c('Blues'=9,'BuGn'=9, 'BuPu'=9, 'GnBu'=9, 'Greens'=9, 'Greys'=9, 'Oranges'=9, 'OrRd'=9, 'PuBu'=9, 'PuBuGn'=9, 'PuRd'=9, 'Purples'=9, 'RdPu'=9, 'Reds'=9, 'YlGn'=9, 'YlGnBu'=9, 'YlOrBr'=9, 'YlOrRd'=9),
    diverging = c('BrBG'=11, 'PiYG'=11, 'PRGn'=11, 'PuOr'=11, 'RdBu'=11, 'RdGy'=11, 'RdYlBu'=11, 'RdYlGn'=11, 'Spectral'=11),
    qualitative = c('Accent'=8, 'Dark2'=8, 'Paired'=12, 'Pastel1'=8, 'Pastel2'=8, 'Set1'=9, 'Set2'=8, 'Set3'=12)
  );

  palettes = unlist(palettes);
  names(palettes) = gsub('\\w+\\.', '', names(palettes))

  if (palette %in% names(palettes))
    i = match(palette, names(palettes))
  else
    i = ((max(c(1, suppressWarnings(as.integer(palette))), na.rm = T)-1) %% length(palettes))+1
    
  col = c();
  col.remain = n;

  while (col.remain > 0)
    {
      if (col.remain > palettes[i])
        {
          next.n = palettes[i]
          col.remain = col.remain-next.n;
        }
      else
        {
          next.n = col.remain
          col.remain = 0;
        }

      col = c(col, brewer.pal(max(next.n, 3), names(palettes[i])))
      i = ((i) %% length(palettes))+1
    }

  col = col[1:n]
  return(col)
}

##################################
# charToDec
#
# converts character vector to byte vector in decimal representation
##################################
charToDec = function(c)
  {    
    return(as(charToRaw(c), 'integer'))
  }

###################################
# which.char
#
# finds the index of the character in subject (length 1 character vector) matching
# nchar = 1 single character query
# eg which.char('a', 'cat') = 2
#
# if query has more than one char (or has length>1) then will return indices matching <any one> of the characters in any
# element of query 
###################################
which.char = function(subject, query)
  {
    if (length(query)>1)
      query = paste(query, collapse = '')
    
    which(charToRaw(subject[1]) %in% charToRaw(query))
  }

##################################
# vaggregate
#
# same as aggregate except returns named vector
# with names as first column of output and values as second
##################################
vaggregate = function(...)
  {
    out = aggregate(...);
    return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
  }


####################################
# modix
#
# converts ix into "modulo index" into a vec of length l
# useful for reusing a vectors contents
####################################
modix = function(ix, l)
  {
    return(((ix-1) %% l)+1)
  }


###############################
# elcycles
#
# enumerates all elementary cycles in a graph via igraph library
#
# A is either an adjacency matrix or igraph object
#
# returns a list 
# $cycles = list of vertices in elementary cycles
# $cycles.eix = list of edges in elementary cycles, where edges are numbered according to the 1D index of adj matrix A
###############################
elcycles = function(A)
  {
    if (inherits(A, 'igraph'))
      A = as.matrix(graph.adjacency(A));
    
    A = abs(sign(A)) * matrix(1:length(A), nrow = nrow(A))
    
    # list of cycles (ie lists of node indices, seeded with self cycles
    out = lapply(which(diag(A)!=0), function(x) x)
    cl = clusters(graph.adjacency(A!=0), 'strong')

    while (length(cl.left <- which(cl$csize>1))>0)
      {                
        # get all cycles
        seeds = match(cl.left, cl$membership)        
        cycles.eix = cycles = list();
               
        for (s in seeds)
          {
            browser()
            tail = which(A[s, ]!=0); # this is vector containing tail item of these.cycles
            these.cycles = lapply(tail, function(x) x)
            these.cycles.eix = lapply(1:length(tail), function(x) list())
            visited = array(FALSE, dim = c(length(these.cycles), ncol(A)))   ## visited is cycles x nodes matrix keeps track of visited nodes in each path
            done = rep(FALSE, length(these.cycles)) ## cycles are done if their tail node comes back to "s"
            while (any(!done)) ## cycle through !done cycles
              {
                print(these.cycles)
                j = which(!done)[1]; # pick next !done path
                i = tail[j]  # i is the last element of that path
                visited[j, i] = TRUE;
                
                if (i != s)
                  {
                    children = which(A[i, ] != 0 & !visited[j, ])
                    children.eix = A[i, children]
                    
                    if (length(children)>0)
                      {
                        these.cycles = c(these.cycles[-j], lapply(children, function(x) c(these.cycles[[j]], x)))
                        these.cycles.eix = c(these.cycles.eix[-j], lapply(children.eix, function(x) c(these.cycles.eix[[j]], x)))
                        visited = rbind(visited[-j, , drop = FALSE], visited[rep(j, length(children)), , drop = FALSE])
                        tail = c(tail[-j], children);
                        done = c(done[-j], done[rep(j, length(children))]);                
                      }
                    else ## remove this path since it is childless and does not end in s
                      {
                        these.cycles = these.cycles[-j]
                        these.cycles.eix = these.cycles.eix[-j]
                        visited = visited[-j, , drop = FALSE]
                        tail = tail[-j];
                        done = done[-j];                        
                      }
                  }
                else
                  done[j] = TRUE ## yes we have found a cycle
              }

            # collect cycles
            cycles = c(cycles, these.cycles);
            cycles.eix = c(cycles.eix, these.cycles.eix)

            # recompute strongly connected components (zeroing out seeds in A matrix);
            A[seeds, ] = 0
            A[, seeds] = 0
            cl = clusters(graph.adjacency(A), 'strong')            
          }        
      }
  }



#########################
# rmix
#
# sample N points from a mixture of k densities of a single functional form (eg norm, beta, multinomial)
# where n is either an integer vector of length k denoting how many samples to be drawn from each density
# (in which case N = sum(n)) or n is a scalar, in which case n points are drawn from each density and N = n*k.
# 
# p = params data frame whose named columns correspond to arguments to rdens (eg $n, $shape1, $shape2 for rbeta or $n, $mean, $sd for rnorm)
# rdens = function encoding random number generator for given density, that takes as input named columns of params
# n = either an nrow(p) integer vector or scalar denoting how many samples to draw from each density 
#
# n can also be just be a column of p
#
# useful for plotting "smears" of points
#
# Output is the rbind-ed output of individual rdens calls
########################
rmix = function(p, rdens, n = NULL)
  {
    if (!is.null(n))
      p$n = n

    tmp = do.call('mapply', c(list(FUN = rdens, SIMPLIFY = FALSE), as.list(p)))
    
    if (length(tmp)>0)
      if (!is.null(dim(tmp[[1]])))
        return(do.call('rbind', tmp))    
      else
        return(do.call('c', tmp))
    else 
      return(tmp)
  }


  
################################
# dmix
#
# generates data frame of density points in a provided range for a provided mix of densities of a given family
# useful for plugging into downstream plotting (eg ggplot 2)
#
# dens is handle of density function
# xlim is a length 2 vector specifying plot bounds 
# n is integer number of points to create per distro
#
# "..." variables depend on density function, arguments should be provided as they would to the
# corresponding R function (ie with respect to vectorization)
#
# dnorm: mean, sd
# dbinom: size, prob
# dmultinom: size, prob
# dgamma: shape, rate
# dbeta: shape1, sheape 2
#
#
# if collapse = TRUE then the densities will be summed according to the mixing parameter yielding a single density
# (ie a fuzzy histogram) summarizing the mixing distribution
################################
dmix = function(dens = 'dnorm', xlim = NULL, n = 500, alpha = NULL, plot = F, fill = T, collapse = F,  ...)
  {
    if (is.null(xlim))
      {
        if (dens == 'dnorm')
          xlim = c(min(list(...)$mean-3*list(...)$sd, na.rm = T), max(list(...)$mean+3*list(...)$sd, na.rm = T))
        else if (dens == 'dbinom')
          xlim = c(0, max(list(...)$size))
        else if (dens == "dbeta")
          xlim = c(0, 1)
        else if (dens == 'dmultinom')
          xlim = c(0, max(list(...)$size))
        else if (dens == "gamma")
          {
            mean = list(...)$alpha/list(...)$beta
            sd = sqrt(list(...)$alpha/list(...)$beta^2)
            xlim = c(mean-3*sd, mean+3*sd)
          }
      }
    
    args = do.call('data.frame', list(...));
    
    if (is.null(alpha))
      alpha = rep(1, nrow(args))/nrow(args)

    x = seq(xlim[1], xlim[2], length = n);
        
    out = data.frame(id = rep(1:nrow(args), each = length(x)),  x = rep(x, nrow(args)), prob = unlist(lapply(1:nrow(args), function(i)
         d = alpha[i]*do.call(dens, c(list(x=x), as.list(args[i, ]))))))

    if (collapse)
      {
        out = aggregate(prob ~ x, data = out, FUN = sum)
        out$id = 1;
      }    
    
    if (plot)
      {
        require(ggplot2)
        if (fill)
          return(ggplot(out, aes(x = x, y = prob, group = id, fill = id)) + geom_ribbon(alpha = 0.3, color = 'black', aes(ymin = 0, ymax = prob)) + scale_fill_gradient(low = 'red'))
        else
          return(ggplot(out, aes(x = x, y = prob, group = id, color = id)) + geom_line(alpha = 0.3) + scale_color_gradient(low = 'red'))
      }
    else  
      return(out)    
  }

#################################
# svec
#
# makes "sparsely" defined numeric vector of length n
# using name= value arguments
#
# eg svec(10, '5'=c(1,2,3,4)) makes vector of length 10 with 1,2,3,4 having value 5.
# Note: numeric keys have to be enclosed in quotes
#
# conflicts (ie values hitting the same index) are resolved with FUN (eg sum = adds)
# (similar to matlab accumarray)
#
# dval is "default" val for unspecified indices
#################################
svec = function(n = 0, dval = 0, op = '+', ...)
  {
    args = list(...);
    n = pmax(n, sapply(args, function(x) max(x, na.rm = T)))
    out = blank = rep(dval, n);
    for (a in names(args))
      {
        tmp = blank;
        tmp[args[[a]]] = as.numeric(a);
        out = do.call(op, list(out, tmp))
      }
    return(out)
  }

##################################
# nz
#
# outputs the nonzero entries of a vector array
##################################
nz = function(x, zero = 0)
  {
    if (length(x)==0)
      return(data.frame())

    if (is.null(nrow(x)))
      x = as.matrix(x)
    
    ix = which(x!=zero, arr.ind = T)

    if (nrow(ix)==0)
      return(data.frame())

    out = NULL;
    
    if (!is.null(rownames(x)))
      out = data.frame(rowname = rownames(x)[ix[,1]], stringsAsFactors = F)

    if (is.null(out))
      out = data.frame(row = ix[,1])
    else
      out = cbind(out, data.frame(row = ix[,1]))
    
    if (!is.null(colnames(x)))
      out = cbind(out, data.frame(colname = colnames(x)[ix[,2]], stringsAsFactors = F))

    out = cbind(out, data.frame(col = ix[,2], val = x[ix], stringsAsFactors = F))
    rownames(out) = NULL
        
    return(out)
  }

 
##################################
# convex.basis
#
# Outputs a matrix K of the convex basis of matrix A
#
# ie each column x = K[,i] is a minimal solution (with respect to sparsity) to 
# Ax = 0, x>=0
#
#
##################################
convex.basis.old = function(A, interval = 80, chunksize = 100)
  {
    ZERO = 1e-8;
    remaining = 1:nrow(A);
    iter = 0;
    i = 0;
#    order = c()
    numelmos = c()
    K_i = I = as(diag(rep(1, ncol(A))), 'sparseMatrix');
#    A_i = as(A %*% K_i, 'sparseMatrix');
    K_i = I = diag(rep(1, ncol(A)))
    A_i = A %*% K_i


    # vector to help rescale matrix (avoid numerical issues)
    mp  = apply(abs(A), 1, min); # minimum value of each column
    mp[mp[ZERO]] = 1; # columns with zero minimum get scale "1"

    st = Sys.time()
    # iterate through rows of A, "canceling" them out
    while (length(remaining)>0)
      {
        iter = iter+1;
        K_last = K_i;

        print(Sys.time() - st)
        
        cat('Iter ', iter, ' Num basis vectors: ', nrow(K_i), " Num active components: ", sum(rowSums(K_i!=0)), "\n")

        i = remaining[which.min(rowSums(A_i[remaining,, drop = FALSE]>=ZERO)*rowSums(A_i[remaining,, drop = FALSE]<=(-ZERO)))]  # chose "cheapest" rows

        remaining = setdiff(remaining, i);
#        order = c(order, i);

        zero_elements = which(abs(A_i[i, ]) <= ZERO);
        K_i1 = K_last[zero_elements, , drop = FALSE];  ## K_i1 = rows of K_last that are already orthogonal to row i of A
        K_i2 = NULL; ## K_i1 = will store positive combs of K_last rows that are orthogonal to row i of A (will compute these below)

        pos_elements = which(A_i[i, ]>ZERO)
        neg_elements = which(A_i[i, ]<(-ZERO))

        cat('Iter ', iter, " Row ", i, ":", length(zero_elements), " zero elements ", length(pos_elements), " pos elements ", length(neg_elements), " neg elements \n")
        
        if (length(pos_elements)>0 & length(neg_elements)>0)
          for (m in seq(1, length(pos_elements), interval))
            for (l in seq(1, length(neg_elements), interval))
              {
                ind_pos = c(m:min(c(m+interval, length(pos_elements))))
                ind_neg = c(l:min(c(l+interval, length(neg_elements))))

                indpairs = cbind(rep(pos_elements[ind_pos], length(ind_neg)),
                  rep(neg_elements[ind_neg], each = length(ind_pos))); # cartesian product of ind_pos and ind_neg
                pix = rep(1:nrow(indpairs), 2)
                ix = c(indpairs[,1], indpairs[,2])
                coeff = c(-A_i[i, indpairs[,2]], A_i[i, indpairs[,1]])
                combs = sparseMatrix(pix, ix, x = coeff, dims = c(nrow(indpairs), nrow(K_last)))
                combs[cbind(pix, ix)] = coeff;

                H = combs %*% K_last;

                # remove duplicated rows in H (with respect to sparsity)
                H = H[!duplicated(as.matrix(H)>ZERO), ];

                # remove rows in H that have subsets in H (with respect to sparsity) .. 
                keep = which(colSums(sparse_subset(abs(H)>ZERO, abs(H)>ZERO, chunksize = chunksize))<=1) # <=1 since every H is its own subset
                H = H[keep, , drop = FALSE]
                
                # remove rows in H that have subsets in K_i2
                if (!is.null(K_i2))
                  {
                    keep = which(colSums(sparse_subset(abs(K_i2)>ZERO, abs(H)>ZERO, chunksize = chunksize))==0) 
                    H = H[keep, , drop = FALSE]
                  }
                
                # remove rows in H that have subsets in K_i1
                if (!is.null(K_i1))
                  {
                    keep = which(colSums(sparse_subset(abs(K_i1)>ZERO, abs(H)>ZERO, chunksize = chunksize))==0) 
                    H = H[keep, , drop = FALSE]
                  }
                
                # maintain numerical stability
                if ((iter %% 10)==0)
                  H = diag(1/apply(abs(H), 1, max)) %*% H
                
#                K_i2 = rBind(K_i2, H)
                K_i2 = rbind(K_i2, as.matrix(H))
              }

#        K_i = rBind(K_i1, K_i2)
        K_i = rbind(K_i1, K_i2)
        A_i = A %*% t(K_i)
        
        if (nrow(K_i)==0)
          return(matrix())
      }
    
    return(t(K_i))       
  }

###############################################
# sparse subset
#
# given k1 x n matrix A and k2 x n matrix B
# returns k1 x k2 matrix C whose entries ij = 1 if the set of nonzero components of row i of A is
# a (+/- strict) subset of the nonzero components of row j of B
#
###############################################
sparse_subset = function(A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
  {
    nz = colSums(A!=0, 1)>0
    if (is.null(dim(A)) | is.null(dim(B)))
      return(NULL)

    C = sparseMatrix(i = c(), j = c(), dims = c(nrow(A), nrow(B)))

    for (i in seq(1, nrow(A), chunksize))
      {
        ixA = i:min(nrow(A), i+chunksize-1)
        for (j in seq(1, nrow(B), chunksize))
          {
            ixB = j:min(nrow(B), j+chunksize-1)

            if (length(ixA)>0 & length(ixB)>0 & !quiet)
              cat(sprintf('\t interval A %s to %s (%d) \t interval B %d to %d (%d)\n', ixA[1], ixA[length(ixA)], nrow(A), ixB[1], ixB[length(ixB)], nrow(B)))
            if (strict)
              C[ixA, ixB] = (sign((A[ixA, , drop = FALSE]!=0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))) * (sign((A[ixA, , drop = FALSE]==0)) %*% sign(t(B[ixB, , drop = FALSE]!=0))>0)
            else
              C[ixA, ixB] = (sign(A[ixA, nz, drop = FALSE]!=0) %*% sign(t(B[ixB, nz, drop = FALSE]==0)))==0
          }
      }
    
    return(C)
  }


######################################################
# morder
#
# matrix order wrt columns ..
# ie ordering rows matrix based on left to right ordering of columns (if MARGIN = 1)
# OR  ordering columns of  matrix based on top to bottom ordering of rows (if MARGIN = 2) 
######################################################
morder = function(A, orient = 1)
  {
    if (orient==1)
      return(do.call('order', lapply(1:ncol(A), function(x) A[,x])))
    else
      return(do.call('order', lapply(1:nrow(A), function(x) A[x,])))
  }


######################################################
# mmatch
#
# match rows of matrix A to matrix B
######################################################
mmatch = function(A, B, dir = 1)
  {
    SEP = ' ';
    Atxt = apply(A, dir, function(x) paste(x, collapse = SEP))
    Btxt = apply(B, dir, function(x) paste(x, collapse = SEP))

    return(match(Atxt, Btxt))
  }

##########################################################
# bisort
#
# "bisorts" matrix according to rows and columns (and optionally removes empty rows, ie with no nonzero)
##########################################################
bisort = function(A, drop = F)
  {
    if (drop)
      A = A[which(rowSums(A!=0)>0), , drop = FALSE]

    A = A[, morder(t(A)), drop = FALSE];
    A = A[morder(A), , drop = FALSE];
    return(A)    
  }

##############################################################
# setxor
#
##############################################################
setxor = function(A, B)  
  {
    return(setdiff(union(A,B), intersect(A,B)))
  }


##############################################################
# sub2ind
#
##############################################################
sub2ind = function(dim, r, c, byrow = F) if (byrow) (r-1)*dim[2] + c else (c-1)*dim[1] + r

##############################################################
# ind2sub
#
##############################################################
ind2sub= function(dim, ind, byrow = F)
  {
    if (byrow)
      cbind(floor((ind-1) / dim[2])+1, ((ind-1) %% dim[2]+1))
    else
      cbind(((ind-1) %% dim[1])+1, floor((ind-1) / dim[1])+1)
  }



#############################################################
# munlist
#
# unlists a list of vectors, matrices, data frames into a n x k matrix
# whose first column specifies the list item index of the entry
# and second column specifies the sublist item index of the entry
# and the remaining columns specifies the value(s) of the vector
# or matrices.
# 
# force.cbind = T will force concatenation via 'cbind'
# force.rbind = T will force concatenation via 'rxsbind'
#############################################################
munlist = function(x, force.rbind = F, force.cbind = F, force.list = F)
  {
    if (!any(c(force.list, force.cbind, force.rbind)))
      {
        if (any(sapply(x, function(y) is.null(dim(y)))))
          force.list = T
        if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1)
          force.rbind = T
        if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1))
          force.cbind = T
      }        
    else
      force.list = T    
        
    if (force.list)
      return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                   iix = unlist(lapply(1:length(x), function(y) if (length(x[[y]])>0) 1:length(x[[y]]) else NULL)), 
                   unlist(x)))
    else if (force.rbind)
      return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))), 
                   iix = unlist(lapply(1:length(x), function(y) if (nrow(x[[y]])>0) 1:nrow(x[[y]]) else NULL)),
                   do.call('rbind', x)))
    else if (force.cbind)
      return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) if (ncol(x[[y]])>0) 1:ncol(x[[y]]) else NULL)),
                   do.call('cbind', x))))
  }


############################################################
# readRDA
# 
# loads Rdata environment into a list variable and returns
# (to mirror RDS functionality)
############################################################
readRDA = function(fn)
  {
    my.env  = new.env()
    load(fn, my.env);
    return(as.list(my.env))
  }


###################################################################
# normalmix.modelsel
#
# model selection for normalmixEM models (in mixtools package)
#
# adapted from repnormmixmodel.sel
#
# k is the maximum model size to entertain
##################################################################
normalmix.modelsel = function(x, k = 2, ...) 
{
  require(mixtools)
  
  if (is.vector(x))
    x = cbind(x)
  
  aic <- NULL
  bic <- NULL
  caic <- NULL
  icl <- NULL
  AIC <- function(emout) {
    emout$loglik - (length(emout$mu) + length(emout$stdev) + 
                    length(emout$lambda) - 1)
  }
  BIC <- function(emout) {
    emout$loglik - log(nrow(x)) * (length(emout$mu) + length(emout$stdev) + 
                                   length(emout$lambda) - 1)/2
  }
  CAIC <- function(emout) {
    emout$loglik - (log(nrow(x)) + 1) * (length(emout$mu) + 
                                         length(emout$stdev) + length(emout$lambda) - 1)/2
  }
  ICL <- function(emout) {
    BIC(emout) - sum(emout$lambda * log(emout$lambda))
  }
  for (i in 1:k) {
    if (i == 1) {
      avx <- as.vector(x)
      mu <- mean(avx)
      s <- sd(avx)
      loglik <- sum(dnorm(avx, mean = mu, sd = s, log = TRUE))
      emout <- list(mu = mu, stdev = s, lambda = 1, loglik = loglik)
    }
    else emout <- normalmixEM(x, k = i, ...)
    aic[i] <- AIC(emout)
    bic[i] <- BIC(emout)
    caic[i] <- CAIC(emout)
    icl[i] <- ICL(emout)
  }
  out = rbind(aic, bic, caic, icl)
  Winner = apply(out, 1, function(x) (1:length(x))[x == max(x)])
  rownames(out) = c("AIC", "BIC", "CAIC", "ICL")
  colnames(out) = 1:k
  cbind(out, Winner)
}

###############################################
# dotplot
#
# Produces dotplot of grouped data 
#
###############################################
dotplot = function(y, group, ylab = '', xlab = '', log = F, dotsize = NULL, binwidth = NULL, title = NULL, ylim = NULL, text.size = NULL)  
  {
    require(ggplot2)
    
    df = data.frame(y = y, group = as.character(group), stringsAsFactors = F)

    if (is.null(dotsize))
      g = ggplot(df, aes(x = group, y = y)) + theme_bw() + theme(text = element_text(size = text.size)) + geom_dotplot(binaxis = 'y', stackdir = 'center', position = 'dodge', binwidth = binwidth)
    else
      g = ggplot(df, aes(x = group, y = y)) + theme_bw() + theme(text = element_text(size = text.size)) + geom_dotplot(binaxis = 'y', stackdir = 'center', position = 'dodge', dotsize = dotsize, binwidth = binwidth)
              
    if (!is.null(ylab))
      g = g + labs(y = ylab)

    if (!is.null(xlab))
      g = g + labs(x = xlab)

    if (!is.null(title))
      g = g + ggtitle(title)
    
    if (log)
      {
        if (!is.null(ylim))
          g = g + scale_y_log10(limits = ylim)
        else
          g = g + scale_y_log10()
      }
    else if (!is.null(ylim))
      g = g + scale_continuous(limits = ylim)                

    return(g)
  }
    


###############################################
# dirr
#
# a variant of dir that gsubs pattern from normal output of dir to name output vector
#
# eg dirr(path, '.txt') will return dir output with .txt removed
# eg dirr(path, '.txt', '.rds' ) will return dir output with .txt subbed with .rds
#
###############################################
dirr = function(x, pattern, rep = '', ...)
  {
    out = dir(x, pattern, ...)
    names(out) = gsub(pattern, rep, dir(x, pattern))
    return(out)
  }
  

###############################
# scrape.http.dir
#
# scrapes links from http directory dump using base.url
# 
###############################
scrape.http.dir = function(base.dir)
  {
    library(RCurl)
    library(XML)
    tmp = as.vector(xpathSApply(htmlParse(getURL(base.dir)), "//a/@href"))
    if (length(tmp)==0)
      return(c())
    else      
      return(paste(base.dir, '/', tmp, sep = ''))
  }


## takes k  n_i x m matrices with n_1, ..., n_k rows and outputs
## a (n_1 * n_2 .. * n_k) x m x k matrix of all k cartesian combinations of
## of the rows of these matrices
##
## first matrix can have 3 dimensions, i.e. be n_1 x m x k_0, in which case
## the additional combinations will be added to the end (i.e. the final
## matrix will havve (n_1 * n_2 * ... * n_k) x  (k_9 + k -1 ) x m combos
.matcart = function(...)
  {
    mats = list(...)
    if (length(mats)==0)
      return(NULL)
    out = mats[[1]]
    if (length(dim(out))==2)
      out = array(out, dim = c(dim(out), 1))
    if (length(mats)==1)
      return(out)
    if (length(mats)>1)
      for (i in 2:length(mats))
        {
          y = mats[[i]]
          ix = cbind(rep(1:nrow(out), nrow(y)), rep(1:nrow(y), each = nrow(out)))
          out = array(c(out[ix[, 1],,], y[ix[,2], ]), dim = c(nrow(ix), dim(out)[2], dim(out)[3]+1))
        }
    return(out)
                                        #        return(aperm(out, c(1, 3, 2)))
  }

###################
# pad 
#
# pads an (integer) vector with k places below and above its lowest and highest value
#
# (by default, clips at 0)
###################
pad = function(x, k, clip = T)
  {
    out = unique(as.vector(rbind(
      apply(cbind(x), 1, function(y) (y-k):(y-1)),
      x,
      apply(cbind(x), 1, function(y) (y+1):(y+k))
      )))

    if (clip)
      out = out[out>0]
    return(out)
  }

##############################
# affine.map
#
# affinely maps 1D points in vector x from interval xlim to interval ylim,
# ie takes points that lie in 
# interval xlim and mapping onto interval ylim using linear / affine map defined by:
# (x0,y0) = c(xlim(1), ylim(1)),
# (x1,y1) = c(xlim(2), ylim(2))
# (using two point formula for line)
# useful for plotting.
#
# if cap.max or cap.min == T then values outside of the range will be capped at min or max
##############################
affine.map = function(x, ylim = c(0,1), xlim = c(min(x), max(x)), cap = F, cap.min = cap, cap.max = cap, clip = T, clip.min = clip, clip.max = clip)
  {
  #  xlim[2] = max(xlim);
  #  ylim[2] = max(ylim);
    
    if (xlim[2]==xlim[1])
      y = rep(mean(ylim), length(x))
    else
      y = (ylim[2]-ylim[1]) / (xlim[2]-xlim[1])*(x-xlim[1]) + ylim[1]

    if (cap.min)
      y[x<min(xlim)] = ylim[which.min(xlim)]
    else if (clip.min)
      y[x<min(xlim)] = NA;
    
    if (cap.max)
      y[x>max(xlim)] = ylim[which.max(xlim)]
    else if (clip.max)
      y[x>max(xlim)] = NA;
    
    return(y)
  }

ppng = function(expr, filename = 'plot.png', height = 1500, width = 1500, cex = 1, ...)
  {
    if (length(cex) == 1)
      cex = rep(cex, 2)
    height = cex[1]*height
    width = cex[2]*width
    filename = paste(normalizePath('~/public_html/'), file.name(filename), sep = '/')
    cat('rendering to', filename, '\n')
    png(filename, height = height, width = width, ...)
    eval(expr)
    dev.off()
  }

ppdf = function(expr, filename = 'plot.pdf', height = 10, width = 10, cex = 1, ...)
  {
    if (length(cex) == 1)
      cex = rep(cex, 2)
    height = cex[1]*height
    width = cex[2]*width
    filename = paste(normalizePath('~/public_html/'), file.name(filename), sep = '/')
    cat('rendering to', filename, '\n')
    pdf(filename, height = height, width = width, ...)
    eval(expr)
    dev.off()
  }


clock = function(expr)
  {
    now = Sys.time()
    eval(expr)
    return(Sys.time()-now)
  }
