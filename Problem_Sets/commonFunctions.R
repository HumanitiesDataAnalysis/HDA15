library(dplyr)
library(ggplot2)
library(tidyr)


folderIntoTD = function (foldername,normalized=F) {
  readSOTU = function(file) {
    message(file)
    text = scan(file,sep="\n",what="raw")
    words = text %>% 
      strsplit("[^A-Za-z]") %>% 
      unlist
    SOTU = data.frame(word=words,filename=file %>% gsub(foldername,"",.) %>% gsub("/","",.) %>% gsub(".txt","",.),stringsAsFactors = FALSE)  
    return(SOTU)
  }
  
  SOTUS = list.files(foldername,full.names = TRUE)
  
  #SOTUS = SOTUS[grep("/19[89]|/20[10]",SOTUS)]
  
  all = lapply(SOTUS,readSOTU)
  allSOTUs = rbind_all(all)
  
  counts = allSOTUs %>% mutate(word=tolower(word)) %>% group_by(filename,word) %>% summarize(count=n()) %>% ungroup
  
  td = counts %>% group_by(word) %>% 
    filter(sum(count)>1000,word!="") %>% 
    ungroup %>% 
    spread(word,count,fill=0)
  
  
  if (normalized) {
    # allow normalization as an argument
    normalized = counts %>% group_by(filename) %>% mutate(ratio = count/sum(count)) %>% ungroup
    norm_data_frame = normalized %>% group_by(word) %>% filter(sum(count)>1000,word!="") %>% select(-count)
    td  = norm_data_frame %>% ungroup %>% spread(key=word,value=ratio,fill=0)
  }
  
  return(td)
  
}


#' Runs a Bookworm query on a local instance
#'
#'@query A bookworm API call, in the form of an R list (will be coerced to JSON)
#'@host the url of the bookworm being queried
#'@method Should just be the default
#'
#'

library(RJSONIO)
bookworm = function(
  host="benschmidt.org",
  database="SOTUgeo2",
  method="return_tsv",
  counttype=list("WordCount"),
  groups = list("year"),
  search_limits = list(),
  query=list()
) {
  for (term in c("method","database","groups","search_limits","counttype")) {
    if(is.null(query[[term]])) {
      query[[term]] = get(term)
    }
    
  }
  if (length(query[['search_limits']]) == 0) {query[['search_limits']]=emptyNamedList}
  if (length(query)==0) {
    query[['database']] = database
  }
  
  require(RCurl)  
  json = toJSON(query)
  json = URLencode(json)
  destination = paste(host,"/cgi-bin/dbbindings.py?query=",json,sep="")
  
  if (method!="return_tsv") {
    data = scan(textConnection(getURL(destination)),what='raw',quote=NULL)
    data = data[data!="===RESULT==="]
    data = paste(data,collapse=" ")
    if (try(assign("data",fromJSON(data[1])),silent=T)==FALSE) {
      warning(destination)
    }
  }
  if (method=="return_tsv") {
    message(destination)
    data = read.table(
      text = getURL(destination),
      header=T,
      sep="\t",
      stringsAsFactors=FALSE,
      blank.lines.skip=T,
      encoding="UTF-8",
      flush=T,
      quote='',
      fill=T,
      comment.char='')
    if(ncol(data)==1 & method=="return_tsv") {
      data = data[grep("^[<>]",data[,1],invert=T),]
      message(destination)
      return(paste(as.character(data),collapse="\n"))
    }
    data[,ncol(data)] = as.numeric(as.character(data[,ncol(data)]))
  }
  return(data)
}






dunning.log = function(set1,set2) {
  wordlist =  merge(set1,set2,by=intersect(names(set1),names(set2)[grep("WordCount",names(set2),invert=T)]),all.x=T,all.y=T)
  #takes a data frame with columns "word," "count.x" and "count.y"
  #Formula (whence variable names) taken from http://wordhoard.northwestern.edu/userman/analysis-comparewords.html
  wordlist$WordCount.x[is.na(wordlist$WordCount.x)] = .5
  wordlist$WordCount.y[is.na(wordlist$WordCount.y)] = .5
  
  wordlist$count.x = wordlist$WordCount.x
  wordlist$count.y = wordlist$WordCount.y  
  attach(wordlist)
  
  wordlist[wordlist==0] = .5
  c = sum(count.x); d = sum(count.y); totalWords = c+d
  wordlist$count.total = count.x+count.y
  wordlist$exp.x = c*(wordlist$count.total)/(totalWords)
  wordlist$exp.y = d*(wordlist$count.total)/(totalWords)
  wordlist$over.x = wordlist$count.x - wordlist$exp.x
  wordlist$over.y = wordlist$count.y - wordlist$exp.y
  
  wordlist$score = 2*(
    (wordlist$count.x*log(
      wordlist$count.x/wordlist$exp.x)) + 
      wordlist$count.y*log(wordlist$count.y/wordlist$exp.y))
  #This makes it negative if the second score is higher
  wordlist$score = wordlist$score * ((wordlist$over.x > 0)*2-1)
  detach(wordlist)
  dunning = wordlist$score
  
  names(dunning) = apply(
    wordlist[,names(set1)[grep("[Cc]ount",names(set1),invert=T),drop=F],drop=F],1,paste,collapse=" ")
  
  data.frame(word = names(dunning),largerIn = ifelse(dunning>0,"set1","set2"),dunning=abs(dunning))
}


