calcFreq <- function(sequence) {
  # sequence = "ACTGTGACTGA...rest of genome"
  out <- list()
  for (i in c("A","C","G","T")) {
    count <- grepl(i,unlist(strsplit(sequence,split="")))
    out[[i]] = sum(count)/length(count)
  }
  return(out)
}

pwm2logodds <- function(pwm,bg.list=list("A"=0.25,"C"=0.25,"G"=0.25,"T"=0.25)) {
  # First make matrix of bg freq
  # listFreq <- list("A"=0.12,"G"=0.12,"C"=0.12,"T"=0.12)
  bg.freq <- matrix(c(bg.list$A,bg.list$C,bg.list$G,bg.list$T),
                    nrow=dim(pwm)[1],ncol=dim(pwm)[2],byrow=T)
  log.odds<-log2(pwm/bg.freq)
  return(log.odds)
}