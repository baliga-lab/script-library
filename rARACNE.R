processARACNEadj <- function(file) {
  # Read in file
  tmp <- readLines(file)
  # Capture parameters
  out <- list()
  for (i in grep(">",tmp)) {
    #print(i)
    l <- strsplit(tmp[i],"   ")[[1]]
    n <- strsplit(l[[1]],"  ")[[1]][[2]]
    n<-sub(" ",".",n)
    n <- gsub("^ +","",n)
    n <- gsub(" +$","",n)
    if (length(l) > 1) {
      if (length(l)==2) {j <- 2} else {j <- 3}
      d <- gsub("^ +","",l[[j]])
      d <- gsub(" +$","",d)
    } else d = ""
    out[[n]] <- d  
  }
  m <- length(tmp)-length(grep(">",tmp))
  cat(paste("Looks like you have", m, "entries. Creating a",m,"-by-",m,"matrix\n"))
  mi.l <- sapply(seq((max(grep(">",tmp))+1),length(tmp)),
                 function(i){i<-strsplit(tmp[i],"\t")})
  mi.names <- sapply(1:length(mi.l),function(i){i<-mi.l[[i]][1]})
  mi.m <- matrix(ncol=m,nrow=m,dimnames=list(mi.names,mi.names))
  diag(mi.m) <- 1
  # Make mask
  mask1<-seq(3,length(mi.l[[1]]),by= 2)
  mask2<-seq(2,length(mi.l[[1]]),by= 2)
  # Parse 
  for (i in seq(1,length(mi.l))) {
    #print(i)
    #entry <- mi.l
    mi.m[mi.l[[i]][1],mi.l[[i]][mask2]] <- as.numeric(mi.l[[i]][mask1])
    }
  out[["Mutual.Information"]] <- mi.mq
  return(out)
}

aracneOnTheFly <- function(dataMatrix,k=0.0907936) {
  # Format data for running with ARACNE
  # ID  Annotation  Cond.i  Cond.j  etc
  # Assuming dataMatrix has col.names and row.names
  cN <- c("ID","Annotation",colnames(dataMatrix))
  rN <- rownames(dataMatrix)
  data <- cbind(rN,rN,dataMatrix)
  colnames(data) <- cN
  # Write tmp file
  write.table(data,file="tmpARACNE.exp",sep="\t",row.names=F,col.names=T,quote=F)
  # Run ARACNE 
  # Default kernel width = 0.0907936
  # Determined from run on full EGRIN2 data set
  # Path to $HOME/ARACNE/usage.txt -- REQUIRED TO RUN
  toL <- gsub("/aracne2","",system("which aracne2",inter=T))
  wd <- getwd()
  #system(paste("ln -s",paste(toL,"/usage.txt",sep=""),paste(wd,"/usage.txt",sep="")))
  # Run ARACNE
  if (!k) {
    system(paste("aracne2 -i tmpARACNE.exp -o tmpARACNEout.adj -H",toL))
  } else {
    # Run with predetermined kernel width
    system(paste("aracne2 -i tmpARACNE.exp -o tmpARACNEout.adj -k",k,"-H",toL))
  }
  out<-processARACNEadj("tmpARACNEout.adj")
  # Clean up
  #system("rm tmpARACNE.exp tmpARACNEout.adj usage.txt")
  return(out)
}
              
              
              