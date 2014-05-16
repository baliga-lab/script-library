# Given an organism name shortcode, this script uses matrix of gene expression ratio
# from Network Portal for a given list of genes from a local file. It uses ggplot2
# library to plot either heatmap, lineplot or smoothed condional mean plot

plotExpression <- function(org=NULL, input=NULL, lineplot=F, heatmap=F, smooth=F){
  require('ggplot2');require('reshape2')
  # get expression ratios file
  ratios.file = as.matrix(read.delim(
    paste("http://networks.systemsbiology.net/static/data-files/expression/",org,"-ratios.tsv", sep=""), header=T, sep="\t"))
  rownames(ratios.file) <- ratios.file[,1]
  # read data file or capture from clipboard
  # we assume inout file has a header and one column
  input.file = read.delim(input, header=T, sep="\t")
  
  # convert to list
  input.genes <- list(as.character(input.file[[1]]))
  #names(input.genes) <- colnames(input.file)
  
  # create ratios data frame for the gene with conditions selected
  tmp1 <- data.frame(ratios.file[unlist(input.genes),], stringsAsFactors=F)
  tmp2 = tmp1[,-c(1,2)]
  exp.cond = dim(tmp2)[2]
  
  # reshape data into longer format
  exp_mat <- melt(t(tmp2))
  exp_mat$value <- as.numeric(as.character(exp_mat$value))
  
  # Lineplot
  if(lineplot){
    p <- ggplot(exp_mat, aes(Var1, value, group=Var2, color=Var2))
    plot <- p + geom_path(stat="identity", position="dodge") +
      theme(axis.text.x = element_blank(), strip.text.x = element_text( angle = 90, size = 8) ) +
      labs(x="Conditions", y="Expression (Log ratios)", title = paste("Expression plot in", exp.cond, "conditions", sep=" " )) +
      guides(col = guide_legend(nrow = 18) ) +
      scale_colour_discrete(name = "Genes")
  }
  
  # Heatmap
  if(heatmap){
    p1 <- ggplot(exp_mat, aes(Var1, Var2, fill=value)) 
    plot <- p1 + geom_tile(aes(fill = value)) +
      scale_fill_gradient(low="blue", high="yellow", space="Lab", name="Log ratios") +
      labs(x="Conditions", y="Genes")
  }
  
  # Smoothed contional mean plot
  if(smooth){
    p2 <- ggplot(exp_mat, aes(Var1, value, group=Var2, color=Var2))
    plot <- p2 + geom_smooth(alpha="0.1") +
      theme(axis.text.x = element_blank(), strip.text.x = element_text( angle = 90, size = 8) ) +
      labs(x="Conditions", y="Smoothed Contional Mean for Expression (Log ratios)", title = paste("Expression plot in", exp.cond, "conditions", sep=" " )) +
      guides(col = guide_legend(nrow = 18)) +
      scale_colour_discrete(name = "Genes")
  }
  return(plot)  
}