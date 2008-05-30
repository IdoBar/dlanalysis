options(showWarnCalls=T, showErrorCalls=T)
if (system("stty -a &>/dev/null") == 0)
  options(width= as.integer(sub(".* ([0-9]+) column.*", "\\1", system("stty -a", intern=T)[1])) - 1 )



util = new.env()

util$unitnorm <- function(x, ...)  (x - mean(x,...)) / sd(x,...)

util$rbern <- function(n, p=0.5)  rbinom(n, size=1, prob=p)

util$msg <- function(...)  cat(..., "\n", file=stderr())

util$strlen <- function(s)  length(strsplit(s,"")[[1]])

util$strmatch <- function(pat,s)  length(grep(pat,s)) > 0

util$strstrip <- function(s)  gsub("^\\s*|\\s*$", "", s)

util$is_empty <- function(collection)  length(collection) == 0

util$as.c <- as.character

util$table.freq <- function(x, ...)  table(x, ...) / sum(table(x, ...))

util$unwhich <- function(indices, len=length(indices)) {
  ret = rep(F,len)
  ret[indices] = T
  ret
}

# util$merge_vec <- function(df, y, by, name) {
#   right = data.frame(bla=y)
#   right[[name]] = right$bla
#   rm(right$bla)
#   right[[by]] = as.numeric(names(y))
#   merge(df, right, sort=FALSE)
# }

util$lax_rbind <- function(...) {
  inputs = list(...)
  each_names = sapply(inputs, names)
  all_names = unique(c(each_names, recursive=TRUE))
  for (k in 1:length(inputs)) {
    if (is.null(inputs[[k]])) next
    more = setdiff(all_names, names(inputs[[k]]))
    inputs[[k]][,more] = NA
  }
  do.call(rbind, inputs)
}

util$fill_bool <- function(bool, true='yes', false='no') {
  ret = rep(NA,length(bool))
  names(ret) = names(bool)
  ret[bool] = true
  ret[!bool] = false
  ret
}

util$trmap <- function(vec, translation_table) {
  ret = rep(NA, length(vec))
  for (x in names(translation_table))
    ret[as.c(vec)==x] = translation_table[[x]]
  ret
}

util$bgrep <- function(pat,x, ...) {
  # "boolean" grep: return a logical vector ready for &, | etc ops.
  # so bgrep works in the world of vector ops like ==, %in%, etc.
  unwhich(grep(pat,x,...), length(x))
}

util$tapply2 <- function(x, ...) {
  # slightly nicer than tapply for some reason, i dont remember
  if (is.factor(x)) {
    r = factor(tapply(as.character(x), ...), levels=levels(x))
  } else {
    r = tapply(x, ...)
  }
  r
}

util$inject <- function(collection, start, fn) {
  acc = start
  for (x in collection)
    acc = fn(acc, x)
  acc
}

util$select <- function(collection, fn) {
  r = c()
  for (x in collection)
    if (fn(x))
      r = c(r, x)
  r
}

util$xprod <- function(xs,ys) {
  ret = list()
  i=0
  for (x in xs)  for (y in ys) {
    i = i+1
    ret[[i]] = list(x=x,y=y)
  }
  ret
}

util$timeit <- function(x) {
  start = Sys.time()
  ret = eval(x)
  finish = Sys.time()
  print(finish - start)
  invisible(ret)
}

util$dotprogress <- function(callback, interval=100) {
  count = 0
  return(function(...) {
    if ((count <<- count+1) %% interval == 0)
      cat(".")
    callback(...)
  })
}

# dataframe-outputting apply and aggregation functions

# like sapply/lapply except it expects fn() to yield lists.
# each list gets coerced into a single row of a dataframe.

util$dfapply <- function(collection, fn, t=TRUE) {
  r = sapply(collection, fn)
  if (t)  r = base::t(r)
  r = matrix2df(r)
  r
}

# sapply() with fn() yielding lists retrns a matrix with named rows/cols ... 
# and whenever you name-index into this thing it return a list ... yuck
# make that shit more normal.

util$matrix2df <- function(x) {
  if (class(x) != 'matrix') stop("why is class ",class(x))
  colnames = names(x[1,])
  data.frame(
    sapply(colnames, function(n) unlist(x[,n])),
    row.names=row.names(x))
}


# like by() but the data types are less crazy:
#  if fn() returns a list, a data frame is returned.  
#    -> byvals are the row names.
#    -> each list is coerced into the rows.
#  if fn() returns a nonlist, a list is returned.
#    -> byvals are the names.
# We attempt to be tolerant for slight inconsistencies in fn()'s return values.

util$dfagg <- function(d, byvals, fn) {
  if (class(byvals) == 'function')
    byvals = byvals(d)
    
  b = by(d, byvals, fn)

  cols = NULL
  for (i in 1:min(100,length(b))) {
    cols = c(cols, names(b[[i]]))
  }
  cols = unique(cols)

  ret = data.frame(row.names=names(b))

  for (col in cols) {
    ret[,col] = sapply(names(b), function(k) b[[k]][[col]])
  }
  if(length(cols) == 0) {
    return(sapply(names(b), function(k) b[[k]]))
  }
  ret
}

util$mymerge <- function(x,y, row.x=F,row.y=F, by=NULL, ...) {
  if (row.x)  x[,by] = row.names(x)
  if (row.y)  y[,by] = row.names(y)

  ret = merge(x,y,by=by, ...)
  if (row.x)  row.names(ret) = row.names(x)
  if (row.y)  row.names(ret) = row.names(y)
  ret
}

# util$flipleft <- function(x,newcol, y,key.y=names(y), key.x=row.names(x)) {
#   right = data.frame(key=key.y)
#   merged = merge(data.frame(key=key.x, x), data.frame(key=key.y,  value=y), by='key')
#   merged$value
#   # merged$value
# }

util$read.xmlss <- function(f) {
  ## BUG: the xml skips cells sometimes.  tricky to parse, argh
  # Mac Excel 2004 calls this "XML Spreadsheet".  It's nice because it's UTF-8.
  #  [ mac .xls seems to be macroman, but xls2csv (perl converter) f's it up,.
  #    and then iconv can't recover.  boo! ]
  csv_pipe = pipe(paste('ruby <<EOF
    require "rubygems"
    require "hpricot"
    require "fastercsv"
    h = Hpricot(File.read("',f,'"))
    mat = (h.at("worksheet")/"row").map{|row| (row/"cell").map{|data| data.inner_text}}
    mat.each{|row| puts row.to_csv}
', sep=''))
  df = read.csv(csv_pipe)
  # close(csv_pipe)
  df
}

########

# for interactivity...

util$excel <- function(d) {
  con = file("/tmp/tmp.csv", "w", encoding="MACROMAN")
  write.csv(d, con)
  system("open -a 'Microsoft Excel' /tmp/tmp.csv")
  close(con)
}

util$mate <- function(...) {
  system(paste("mate", ...))
}

# pretty-print as yaml.  intended for rows with big textual cells.

util$ppy <- function(x, column.major=FALSE, ...) {
  library(yaml)
  cat(as.yaml(x, column.major=column.major), ...)
  cat("\n", ...)
}


while("util" %in% search())
  detach("util")
attach(util)