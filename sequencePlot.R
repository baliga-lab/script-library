nucleotide2color <- function(sequence) {
  # Take nucleotide sequence converts to colors
  if (class(sequence)=="character") {
    seq <- unlist(strsplit(sequence,split=""))
  }
  map <- list("A"="green","a"="green",
              "C"="blue","c"="blue",
              "G"="orange","g"="orange",
              "T"="red","t"="red",
              "N"="grey","n"="grey")
  o <- unlist(map[seq])
  return(o)
}

plotGeneMotifNucleotideSequence <- function(sequence,motifMatrix,trim=T,displayAll=F,heightCutoff=NULL,seqName=NULL...) {
  # Sequece is character of vector of sequence elements
  # motifMatrix is a matrix of pssm values or any value
  # that can be translated to height of sequence
  if (!is.null(seqName)) {
    names(sequence) = seqName
  }
  if (!displayAll) {
    if (!is.null(heightCutoff)) {
      rnames <- rownames(motifMatrix)
      ind <- which(apply(motifMatrix,1,max)>=heightCutoff)
      motifMatrix <- as.matrix(motifMatrix[ind,])
      if (dim(motifMatrix)[2]==1) { #Column vector
        motifMatrix <- t(motifMatrix)
        rownames(motifMatrix) <- rnames[ind]
      }
    } else {
      rnames <- rownames(motifMatrix)
      ind <- which(apply(motifMatrix,1,max)>0)
      motifMatrix <- as.matrix(motifMatrix[ind,])
      if (dim(motifMatrix)[2]==1) { #Column vector
        motifMatrix <- t(motifMatrix)
        rownames(motifMatrix) <- rnames[ind]
      }
    }
  }
  if (trim) {
    start <- min(which(motifMatrix>0,arr.ind=T)[,2])
    sequence <- substr(sequence,start,nchar(sequence))
    motifMatrix <- as.matrix(motifMatrix[,start:dim(motifMatrix)[2]])
    if (dim(motifMatrix)[2]==1) { #Column vector
        motifMatrix <- t(motifMatrix)
        rownames(motifMatrix) <- rnames[ind]
      }
  }
  if (dim(motifMatrix)[1]>15) {
    print("You have too many motifs. Only plotting 15.")
  }
  cex <- par()$cex
  colors <- nucleotide2color(sequence)
  yMax = max(motifMatrix)
  par(mfrow=c(16,1),mar=c(1,4,1,.5))
  barplot2(rep(1,nchar(sequence)),space=0,xlim=c(0,nchar(sequence)),ylim=c(0,5),col=colors,plot.grid=F,axes=F,main=names(sequence),names.arg=F)
  par(mar=c(.5,4,1,.5))
  for (i in 1:dim(motifMatrix)[1]) {
    barplot2(motifMatrix[i,],space=0,xlim=c(0,dim(motifMatrix)[2]),col=colors,plot.grid=T,ylab = rownames(motifMatrix)[i],ylim=c(0,yMax),names.arg=F,cex.axis=cex*.6,cex.lab=cex*0.8)
  }
}