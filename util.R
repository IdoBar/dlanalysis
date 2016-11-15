# util.R:
# Utilities to make R a happier place
# Brendan O'Connor, brenocon.com/code


########################################
## Put everything into an environment, to not pollute global namespace

util = new.env()

########################################
## Better I/O routines
util$repath <- function(x) {
   xa <- gsub('\\\\', '/', x)
 }

util$read.tsv <- function(..., header=F, sep='\t', quote='', comment='', na.strings='', stringsAsFactors=FALSE) {
  # read.table() wrapper with default settings for no-nonsense, pure TSV
  # Typical use case is output from another program.
  # (R's defaults are more geared for human-readable datafiles, which is less
  # feasible for large-scale data anyway.)
  # These options are substantially faster than read.table() defaults.
  #   (see e.g. LINK)
  # stringsAsFactors is the devil.

  args = list(...)
  args$header = header
  if (!is.null(args$col.names)) {
    # read.delim() is not smart about this.  Yikes.
    args$header = FALSE
  }
  args$sep = sep
  args$quote = quote
  args$comment = comment
  args$stringsAsFactors = stringsAsFactors
  args$na.strings = na.strings
  do.call(read.delim, args)
}

util$write.tsv <- function(..., header=NA, col.names=F, row.names=F, sep='\t', na='', quote=F) {
  # 'header' to 'col.names' naming consistency with read.table()
  if (is.finite(header)) col.names = header
  write.table(..., col.names=col.names, row.names=row.names, sep=sep, na=na, quote=quote)
}

# write fasta files
util$write.fasta <- function(id_seq_table, line_width=60, filename){
  file.create(filename)
  apply(id_seq_table, 1, function(entry) cat(sprintf("%s\n%s\n\n", entry[1],  
                                           gsub(sprintf("([GTAC]{%i})",line_width), "\\1\n", entry[2])), 
                                             file = fasta_file, append = TRUE))
}
########################################
##  Misc small routines
  # calculates standard error
util$se <- function(x) sqrt(var(x)/length(x))
util$as.c <- as.character

util$unwhich <- function(indices, len=length(indices)) {
  # reverse of which(): from indices to boolean mask.
  ret = rep(F,len)
  ret[indices] = T
  ret
}

util$nna <- function(...) !is.na(...)   # i type this a lot, i think its worth 3 characters + shift key

util$shuffle <- function(...) UseMethod("shuffle")

util$shuffle.default <- function(x)  x[order(runif(length(x)))]

util$shuffle.data.frame <- function(x)  x[order(runif(nrow(x))),]

util$sample_df <- function(d, size=10, ...)  {
  samp = sample(1:nrow(d), size=size, ...)
  d[samp,]
}

util$present_levels <- function(x) intersect(levels(x), x)

util$trim_levels <- function(...) UseMethod("trim_levels")

util$trim_levels.factor <- function(x)  factor(x, levels=present_levels(x))

util$trim_levels.data.frame <- function(x) {
  for (n in names(x))
    if (is.factor(x[,n]))
      x[,n] = trim_levels(x[,n])
  x
}


# grep() returns indices of matches.  Variants:

util$bgrep <- function(pat,x, ...) {
  # "boolean" grep: return a logical vector ready for vector ops
  # like & |  and others
  unwhich(grep(pat,x,...), length(x))
}

util$ngrep <- function(pat,x, ...)
  # "normal" grep: return values, not indices
  x[grep(pat,x,...)]


########################################
##  Other data manipulation routines
  
  # Shortens and "prettify" long terms 
util$pretty_term <- function(term, comma=TRUE, dot=TRUE, paren=TRUE, brace=TRUE, cap=TRUE, pref="PREDICTED:", charNum=50, flag_change=FALSE) {
  term <- sapply(term, function (s) {
    s_ <- s
    if (!is.null(pref)) s_ <- sub(paste0("^",pref, "\\s*"), "", s_)
    if (cap) s_ <- sub("(^.)", "\\U\\1", s_, perl=TRUE)
    if (paren) s_ <- sub(" \\(.+\\)", "",  s_)
    if (brace) s_ <- sub(" \\[.+\\]", "",  s_)
    if (dot) s_ <- unlist(strsplit(s_, ".", fixed = TRUE))[1]
    if (comma) s_ <- unlist(strsplit(s_, ",", fixed = TRUE))[1]
    # Trim to the end of the word if character number limit is reached and remove any trailing commas or dots
    if (charNum>0) s_ <- sub("[,|\\.]$", "", sub(sprintf("(^.{%s}[\\S]*)[ ]*.*", charNum), "\\1", s_, perl=TRUE), perl = TRUE)
    # Add a ^ at the end if the term was trimmed anyhow
    if (flag_change && nchar(s_)<nchar(s) && !grepl("\\*$", s_)) return(paste0(s_, "^"))
    else return(s_)
  }, USE.NAMES = FALSE)
  return(term)
}
  

util$tapply2 <- function(x, ...) {
  # like tapply but preserves factors
  if (is.factor(x)) {
    r = factor(tapply(as.character(x), ...), levels=levels(x))
  } else {
    r = tapply(x, ...)
  }
  r
}

util$inject <- function(collection, start, fn) {
  # like lisp reduce.  (named after ruby)
  acc = start
  for (x in collection)
    acc = fn(acc, x)
  acc
}




########################################
## Printing, viewing
## see also:  str()

util$printf <- function(...) cat(sprintf(...))

util$listprint <- function(x) {
  s = paste(sapply(names(x), function(n)  sprintf("%s=%s", n,x[[n]])), collapse=' ')
  printf("%s\n", s)
}

util$msg <- function(...)  cat(..., "\n", file=stderr())

util$h = utils::head

util$ppy <- function(x, column.major=FALSE, ...) {
  # pretty-print as yaml.  intended for rows with big textual cells.
  # a la mysql's \G operator
  library(yaml)
  cat(as.yaml(x, column.major=column.major), ...)
  cat("\n", ...)
}

util$table_html = function(...) {
  # Intended for inside dosink()
  columns = list(...)
  ncol = length(columns)
  nrow = length(columns[[1]])
  # assume columns are in parallel
  printf("\n<table cellpadding=3 border=1 cellspacing=0 bordercolor=gray>")
  for (i in 1:nrow) {
    printf("\n<tr>")
    for (j in 1:ncol)
      printf("\n  <td>%s", columns[[j]][i])
  }
  printf("\n</table>\n")
}


########################################
##  Workspace management

  # install and load packages from CRAN, Bioconductor or git, see usage below
  
util$install.deps <- function(p, repo="cran"){
  call_install <- switch(repo,cran=c("install.packages(package, repos=\"http://cloud.r-project.org/\""),
                         bioc=c("biocLite(package, suppressUpdates=TRUE"),
                         git=c("install.github(package"))
  if (repo=="bioc") eval(parse(text = getURL("http://bioconductor.org/biocLite.R", ssl.verifypeer=FALSE)))
  for (package in p) {
    if (!package %in% utils::installed.packages()) {
      cat(sprintf("Please wait, installing and loading missing package %s and its dependencies from %s.\n",
                  package, repo), file=stderr())
      suppressWarnings(eval(parse(text = sprintf("%s, quiet = TRUE)",call_install))))
      if (!package %in% utils::installed.packages()) {
        ifelse(!interactive(),q("no", 2, TRUE),
               stop(sprintf("Unable to install missing package %s from %s.\n",
                            package, repo), call. = FALSE))
      }
    }
    require(package, character.only=T, quietly=TRUE, warn.conflicts=TRUE)
  }
}

# CRAN_packages <- c("RCurl", "dplyr", "tidyr", "ggplot2", "RColorBrewer", "RSQLite")
  # install.deps(CRAN_packages)

  
# improved list of objects
# http://stackoverflow.com/questions/1358003/tricks-to-manage-the-available-memory-in-an-r-session
util$list_objects = function (pos = 1, pattern) {
    napply = function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names = ls(pos = pos, pattern = pattern)
    N = length(names)
    obj_class = napply(names, function(x) as.character(class(x))[1])
    obj_mode = napply(names, mode)
    obj_type = ifelse(is.na(obj_class), obj_mode, obj_class)
    obj_prettysize = napply(names, function(x) {
                           capture.output(print(object.size(x), units = "auto")) })
    obj_size = napply(names, object.size)
    obj_prettysize[obj_size < 1e6] = ""

    obj_length = napply(names, function(x) length(x))
    obj_dim = t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))

    is_flat = is.na(obj_dim)[, 1]
    is_vector = napply(names, function(x) is.vector(x) & class(x) != 'list')


    info_width = max(20, options('width')$width - 60)

    small_str = function(x) {
      out = capture.output(
        str(x, max.level=0, give.attr=F, give.head=F, width=info_width, strict.width='cut')
      )
      out = str_c(out,collapse=' ')
      out = cutoff(str_replace(out,"\n"," "))
      if (str_detect(out, "^List of"))
        out = str_c("[Names] $ ", str_c(names(x),collapse=' '))
      cutoff(out)
    }

    cutoff = function(s) {
      if (str_length(s) >= info_width) {
        str_c(str_sub(s,1,info_width-2),'..')
      } else {
        s
      }
    }
      
    pad = function(s) sprintf(" %s", s)

    out <- data.frame(
      Type = obj_type,
      Size = obj_prettysize,
      Dim = ifelse(is_vector | is_flat, obj_length, 
        sprintf("(%s, %s)", obj_dim[,1], obj_dim[,2])),
      Value = napply(names, function(x) 
        if (class(x) %in% c('data.frame','list') && !is.null(names(x)))
          cutoff(str_c("[Names] $ ",str_c(names(x), collapse=' ')))
        else small_str(x)
        ),
      stringsAsFactors=F)
    row.names(out) = names
    out$Dim = sprintf(" %s", out$Dim)
    out$Value = sprintf(str_c(" %-", info_width, "s"), out$Value)

    out = rbind(subset(out, Type!='function'), subset(out, Type=='function'))
    out
}

util$lsos = function() {
  d = list_objects()
  d$name = row.names(d)
  d = subset(d, name != 'util')
  row.names(d)=d$name
  d$name=NULL
  d
}


########################################
## For performance optimization and long-running jobs

util$timeit <- function(expr, name=NULL) {
  # print how long the expression takes, and return its value too.
  # So you can interpose   timeit({ blabla })   around any chunk of code "blabla".
  start = Sys.time()
  ret = eval(expr)
  finish = Sys.time()
  if (!is.null(name)) cat(name,": ")
  print(finish-start)
  invisible(ret)
}

util$dotprogress <- function(callback, interval=10) {
  # intended to wrap the anonymous callback for sapply() or somesuch.
  # ALTERNATIVE: plyr *ply(.progress='text')
  count = 0
  return(function(...) {
    if ((count <<- count+1) %% interval == 0)
      cat(".")
    callback(...)
  })
}


########################################
##  External programs for interactivity



########################################
## Graphics output wrappers
## For easy one-liners, like:
## dopdf("tmp.pdf",width=5,height=5,cmd=plot(x,y))

util$dopdf <- function(filename,..., cmd) {
  pdf(filename, ...)
  eval(cmd)
  dev.off()
  if (exists('OPEN') && OPEN)
    system(sprintf("open %s", filename))
}

util$dopng <- function(filename,..., cmd) {
  png(filename, ...)
  eval(cmd)
  dev.off()
  if ((exists('OPEN') && OPEN))
    system(sprintf("open %s", filename))
}

util$dosink <- function(filename,cmd, open=NULL) {
  # like capture.output() but follows open/OPEN conventions here
  sink(filename)
  eval(cmd)
  sink(NULL)
  if (prio_check(open, exists('OPEN') && OPEN))
    system(sprintf("open %s", filename))
}

util$dosvg <- function(filename, ..., cmd, open=NULL) {
  library("RSvgDevice")
  devSVG(filename, ...)
  eval(cmd)
  dev.off()
  if (prio_check(open, exists('OPEN') && OPEN))
    system(sprintf("open %s", filename))
}


########################################
## Plotting routines

util$linelight <- function(x,y, lty='dashed', col='lightgray', ...) {
  # highlight a point with lines running to the axes.
  left = par('usr')[1]
  bot = par('usr')[3]
  segments(left,y, x,y, lty=lty, col=col, ...)
  segments(x,bot,  x,y, lty=lty, col=col, ...)
}

util$hintonplot <- function(mat, max_value=max(abs(mat)), mid_value=0, ...) {
  # Plots a matrix as colored, size-varying boxes
  # I dunno who started calling this a "Hinton plot", but anyways

  # Example:
  # hintonplot(matrix(rnorm(100),10))

  # Example, for counts:
  # table(cyl=mtcars$cyl, mpg=cut(mtcars$mpg,3))
  #    mpg
  # cyl (10.4,18.2] (18.2,26.1] (26.1,33.9]
  #   4           0           6           5
  #   6           2           5           0
  #   8          12           2           0
  # hintonplot(table(cyl=mtcars$cyl, mpg=cut(mtcars$mpg,3)))

  plot.new()
  plot.window(xlim=c(0.5,ncol(mat)+0.5), ylim=c(0.5,nrow(mat)+0.5))

  x_mid = 1:ncol(mat)
  y_mid = 1:nrow(mat)

  area = abs(mat) / max_value
  side = sqrt(area)

  for (x in 1:ncol(mat)) {
    for (y in nrow(mat):1) {
      # ym = (nrow(mat):1)[y]
      ym = y
      d = side[ym,x] / 2
      rect(x-d, y-d, x+d, y+d, col=if (mat[ym,x]>0) 'darkblue' else 'darkred')
    }
  }

  axis(1, 1:ncol(mat), labels=colnames(mat))
  # axis(2, nrow(mat):1, labels=row.names(mat))
  axis(2, 1:nrow(mat), labels=row.names(mat))
  title(xlab=names(dimnames(mat))[2], ylab=names(dimnames(mat))[1], ...)
}



########################################



########################################
## Has to be last in file

while("util" %in% search())
  detach("util")
attach(util)

