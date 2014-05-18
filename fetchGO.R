dlf <- function (f, url, msg = NULL, mode = "wb", quiet = F, ...) 
{
  err <- 0
  if (mode == "wb" || !file.exists(f) || file.info(f)$size == 
    118) {
    if (!file.exists(dirname(f))) 
      try(dir.create(dirname(f), recursive = T))
    if (!is.null(msg)) 
      cat(msg, "\n")
    err <- try(download.file(url, destfile = f, mode = mode, 
                             quiet = quiet, ...))
  }
  closeAllConnections()
  err
}


load_go_microbes_online <- function(IdOverride=NULL) {
  require(topGO)
  # Currently requires an active egrin env for the species of interest
  cat("Using GO annotations from MicrobesOnline...\n")
  cat("Storing results in ./data/...\n")
  try(dir.create("./data"))
  if (!is.null(IdOverride)) {
    fname <- paste("data/", e$rsat.species, "/microbesonline_geneontology_", 
                   IdOverride, ".named", sep = "")
    err <- dlf(fname, paste("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", 
                            IdOverride, ";export=tab", sep = ""),mode="wb")
  } else {
    fname <- paste("data/", e$rsat.species, "/microbesonline_geneontology_", 
                   e$taxon.id, ".named", sep = "")
    err <- dlf(fname, paste("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=", 
                            e$taxon.id, ";export=tab", sep = ""),,mode="wb")
    if (e$genome.info$org.id$V1[1] != e$taxon.id && (!file.exists(fname) || file.info(fname)$size == 
      118)) {
      fname <- paste("data/", e$rsat.species, "/microbesonline_geneontology_", 
                     e$genome.info$org.id$V1[1], ".named", sep = "")
      err <- dlf(fname, paste("http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?oId=", 
                              e$genome.info$org.id$V1[1], ";export=tab", sep = ""),mode="wb")
    }
  }
  if (file.exists(fname)) 
    cat("Succesfully fetched GO annotations. Parsing...\n")
  f <- read.delim(fname)
  # try to match appropriate names
  # use accession to pull out names that overlap with ratios matrix
  # remove entries without accession
  f <- f[which(sapply(f[,"accession"],nchar)>1),]
  syns <- e$get.synonyms(f[,"accession"])
  syns.trans <- lapply(seq(1,length(syns)),function(i){syns[[i]][syns[[i]]%in%rownames(e$ratios[[1]])][1]})
  ind <- which(sapply(syns.trans,length)>0)
  fname.map <- paste("data/", e$rsat.species, "/microbesonline_geneontology_", 
                 e$genome.info$org.id$V1[1], ".map", sep = "")
  write.table(data.frame(unlist(syns.trans[ind]),f[ind,"GO"]),fname.map,sep="\t",quote=F,row.names=F,col.names=F)
  gene2go <- readMappings(fname.map)
  return(gene2go)
}

load_topgo_map <- function(file) {
  require(topGO)
  gene2go <- readMappings(file)
  return(gene2go)
}

get_topGO_object <- function(genes,gene2go,ontology=c("BP","MF","CC")[1]) {
  require(topGO)
  # genes is a vector containing genes of interest
  geneList <- factor(as.integer(names(gene2go)%in%genes))
  names(geneList) <- names(gene2go)
  GOdata <- new("topGOdata", ontology = ontology, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene2go)
  #GOdata can be used directly for analysis, e.g. 
  # test <- runTest(GOdata,algorithm="classic",statistic="fisher")
  # results <- GenTable(GOdata,test,topNodes=10)
  return(GOdata)
}
