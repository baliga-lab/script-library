writeFASTA <- function(getsequencesList,file) {
  # Needs to be written as a more generic utility
  # For the moment, it requires output from
  # function getsequences()
  # Which is built into the e environment 
  if (file.exists(file)) {
    print(paste("WARNING: File exists. Overwriting",file))
    file.remove(file)
  }
  out<-lapply(1:length(getsequencesList),function(i){
    at<-attr(g.seq,"start.stops")[i,]
    write(paste(
      paste(">",names(getsequencesList)[i],at[[1]],at[[2]],at[[3]],at[[4]],sep=" "),
      "\n",getsequencesList[[i]]),file=file,append=TRUE)
  })
  invisible(out)
}