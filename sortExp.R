sortExp <- function(ratMat,index) {
  # Sort expression of ratMat based on expression of index
  if (index=="mean") {
   m <- colMeans(ratMat)
   o <- order(m)
  } else {
    o <- order(ratMat[index,])
  }
  out<-ratMat[,o]
  return(out)
}