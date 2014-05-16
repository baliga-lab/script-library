# Given an input file with a list of genes this script will calculate
# enrichment of network modules from Network Portal by using hypergeometric
# p-values

geneSetEnrichment <- function(org=NULL, input=NULL, output.file=F){
  
  # grab organism module-to-genes file
  module.to.genes <- read.delim(
    paste("http://networks.systemsbiology.net/", org, "/modgenes/export", sep=""), header=T, sep="\t")
  
  # collect list of genes for each module
  module.members = list()
  module.members = sapply(module.to.genes[,2], function(i){
    strsplit(as.character(i), split=":")
  })    
  names(module.members) <- as.character(module.to.genes[,1])
  
  # get the total number of genes in all modules
  total.genes = length(unique(rle(unlist(module.members))$values))
  
  # read data file or capture from clipboard
  # we assume inout file has a header and one column
  input.file = read.delim(input, header=T, sep="\t")
  
  # convert to list
  input.genes <- list(as.character(input.file[[1]]))
  names(input.genes) <- colnames(input.file)
  
  # function to get hypergeometric pValues
  hyper.pvalues <- function(){
    pvalues = data.frame(stringsAsFactors=F)
    for(input.j in names(input.genes)) {
      for(module.j in names(module.members)) {
        pvalues = rbind(pvalues, cbind(Input.Name = input.j,
                                       Enriched.Module=module.j,
                                       Number.Input.Genes= length(input.genes[[input.j]]),
                                       Number.Module.Genes= length(module.members[[module.j]]),
                                       overlap=length(intersect(input.genes[[input.j]], module.members[[module.j]])),
                                       overlap.genes = paste(intersect(input.genes[[input.j]], module.members[[module.j]]), sep="", collapse=":"),
                                       Enrichment.p.value=phyper(length(intersect(module.members[[module.j]],input.genes[[input.j]])), # q
                                                                 length(input.genes[[input.j]]), # m
                                                                 total.genes - length(input.genes[[input.j]]), # n
                                                                 length(module.members[[module.j]]), # k
                                                                 lower.tail=F) ) )
        cat("Input:", input.j, "vs ", "Module:", module.j, "\n" )
      }
    }
    return(pvalues)
  }
  pValues <- hyper.pvalues()
    
  # filter hypergeomtric pvalues and correct for multiple comparison
  filter.pvalues <- function(){
    ff.pvalues = cbind(pValues, Corrected.p.value = p.adjust(as.numeric(as.matrix(pValues)[,'Enrichment.p.value']), method='BH'))
    fil1.pvalues = ff.pvalues[which(ff.pvalues$Corrected.p.value < 0.005),]
    cat("There are", dim(fil1.pvalues)[[1]], "comparison with BH p.value smaller than 0.005 \n")
    return(fil1.pvalues)
  }
  filtered.pvalues <- filter.pvalues()
  return(filtered.pvalues)
  # write into file if output.file is TRUE
  if(output.file){
    write.table(filtered.pvalues, file="gsea_result.txt", row.names=F, sep="\t")
    cat("Results were written into file gsea_result.txt \n")
  }
}