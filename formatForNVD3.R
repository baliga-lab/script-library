formatForNVD3 <- function(data_list,type=c("bar","scatter"),...){
		# compile all vs your score for a measure
		# format for JSON and plotting in nvd3
		# data_list is a list (ordered) of data frames containing (at a minimum) x and y values
		all = cbind(all$breaks,all$density); colnames(all) <- c("x","y")
		all = data.frame(all)
		yours <- list(data.frame(x=data_scaled[length(data_scaled)],y=max(all$y)))
		x <- data.frame(key = c("All", "User"))
		x$values = list(all,yours)
		return(x)
}