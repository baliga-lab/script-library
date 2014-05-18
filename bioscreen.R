####################################################################################
# Author: Aaron Brooks
# Affiliation: Institute for Systems Biology, Seattle, WA
# Date of creation: 03/03/2014
# Last update: 03/03/2014
####################################################################################
# DESCRIPTION
####################################################################################
# functions for parsing bioscreen data
#################################################################################### 

loadBioscreen <- function(file,sep=",",time_format = c("min","hours")[1],exclude){
	data = read.delim(file,sep,header=T)
	# parse time to minutes
	parseTime <- function(time,format="min") {
		sapply(time,function(x) {
			y = as.numeric(strsplit(as.character(x),split=":")[[1]])
			if (format=="min") {
				return(y[1]*60 + y[2] + y[3]/60)
			} else if (format=="hours") {
				return(y[1] + y[2]/60 + y[3]/(60*60))
				} else {
					return(x)
				}
			})
	}
	parseColumnNames <- function(n) {
		# return name:indexes
		# find unique names - ignore replicates
		uniqueN <- lapply(n,function(x) {
			nsplit = strsplit(x,"\\.")[[1]]
			if (length(nsplit)>1) {
				if (is.na(as.numeric(nsplit[length(nsplit)]))==T) {
					nsplit = paste(nsplit[1:length(nsplit)],collapse=".")
				} else {
					nsplit = paste(nsplit[1:length(nsplit)-1],collapse=".")
				}
			} 
			return(strsplit(nsplit,split="_")[[1]])
			})
		uniqueN_list <- lapply(unique(unlist(uniqueN)),function(x){
				unlist(sapply(seq(1,length(uniqueN)),function(j){
					if (length(grep(paste("^",x,"$",sep=""),uniqueN[[j]]))>0) {
						return(j)
					}
				}))
			}); names(uniqueN_list) <- unique(unlist(uniqueN))
	
		return(uniqueN_list)
	}
	QC <- function(data,nameList) {
		# quick QC to identify suspicious data

	}
	time = parseTime(data[,1],time_format)
	data = data[,-1]
	if (length(exclude)>0) {
		for (i in exclude) {
			data = data[,-i]
		}
	}
	nameList = parseColumnNames(colnames(data))
	to.r <- list()
	to.r$time = time
	to.r$data = data; rownames(to.r$data) = time
	to.r$replicate.mapping = nameList
	return(to.r)
}

combine <- function(pattern=c("M9","Kan","Arg"),mapping) {
	to.r <- sort(unique(unlist(data$replicate.mapping)))
	for (i in pattern) {
		to.r <- intersect(to.r,mapping[[i]])
	}
	return(to.r)
}

selectMultiple <- function(indexes,data,dir = c(1,2)[2]) {
	if (dir==1) {
			to.r <- do.call(cbind,lapply(indexes,function(i){
				data[i,]
				})) 
			colnames(to.r) <- colnames(data)
			rownames(to.r) <- rownames(data)[indexes]
		} else if (dir==2) {
			to.r <- do.call(cbind,lapply(indexes,function(i){
				data[,i]
				})) 
			colnames(to.r) <- colnames(data)[indexes]
			rownames(to.r) <- rownames(data)
		}
	return(to.r)
}

