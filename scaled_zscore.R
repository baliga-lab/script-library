scaled_zscore <- function(input) {
  if (is.null(dim(input))) {
    input.sd <- sd(input,na.rm=T)
  } else {
    input.sd<-sd(as.vector(input),na.rm=T)
  }
  input.mean <- mean(input,na.rm=T)
  zscore <- (input-input.mean)/input.sd
  max.z <- max(zscore,na.rm=T)
  min.z <- min(zscore,na.rm=T)
  scaled.zscore <- (zscore-min.z)/(max.z-min.z)
  return(scaled.zscore)
}

zscore <- function(input) {
  if (is.null(dim(input))) {
    input.sd <- sd(input,na.rm=T)
  } else {
    input.sd<-sd(as.vector(input),na.rm=T)
  }
  input.mean <- mean(input,na.rm=T)
  zscore <- (input-input.mean)/input.sd
  return(zscore)
}

scaled <- function(input) {
  max.input <- max(input,na.rm=T)
  min.input <- min(input,na.rm=T)
  scaled.input <- (input-min.input)/(max.input-min.input)
  return(scaled.input)
}