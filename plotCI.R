plotCI <- function(data,upperCI,lowerCI=NULL,cols=NULL,sort=F,autoY=T,autoX=T,stack=F,legend=T,position="bottomright",...) {
  ####################################################################################
  # Author: Aaron Brooks
  # Affiliation: Institute for Systems Biology, Seattle, WA
  # Date of creation: 11/10/2011
  # Last update: 03/12/2012
  ####################################################################################
  # DESCRIPTION
  ####################################################################################
  # This function plots data with confidence bounds using polygons
  # Useful for visualizing trends in groups of data without the
  # noise of many lines
  # DATA is a vector or list of data to be plotted
  # The names() attribute of each element in DATA specifies 
  # x coordinates (if applicable)
  # UPPERCI is a vector or list of distances to be plotted 
  # upward from data
  # same goes for LOWERCI, but this is plotted on bottom
  # lengths of DATA, UPPERCI, LOWERCI should correspond
  # COLS allows for custom colors
  # SORT controls sorting of data
  # AUTOLIM will automatically set limits if TRUE
  # STACK will plot list results side-by-side rather than on top 
  # of one another
  ####################################################################################
  require(RColorBrewer)
  if (F) {
    # EXAMPLES
    # vector: 
    data=as.numeric(sample.int(100));upperCI=as.numeric(sample.int(100)/10);lowerCI=as.numeric(sample.int(100)/10)
    plotCI(data,upperCI,lowerCI,sort=F,xlab="Sorted Index",ylab="Value of Random Number",main="Random Vector Input, Unsorted")
    plotCI(data,upperCI,lowerCI,sort=T,xlab="Sorted Index",ylab="Value of Random Number",main="Random Vector Input, Sorted")
    # vector time series:
    # noisy sine wave
    data = sapply(seq(1:100),sin)
    names(data) = seq(1:100)
    upperCI = sample(seq(0,1,length.out=1000),100)
    lowerCI = sample(seq(0,1,length.out=1000),100)
    plotCI(data,upperCI,lowerCI,sort=F,xlab="x",ylab="sin(x)",main="Time Series Data")
    # list: 
    data=lapply(seq(1:3),function(i){sample.int(100)})
    upperCI=lapply(seq(1:3),function(i){sample.int(100)/10})
    lowerCI=lapply(seq(1:3),function(i){sample.int(100)/10})
    plotCI(data,upperCI,lowerCI,sort=T,xlab="Sorted Index",ylab="Value of Random Number",main="List Input")
    # list stacked:
    data=lapply(seq(1:3),function(i){sample.int(100)})
    upperCI=lapply(seq(1:3),function(i){sample.int(100)/10})
    lowerCI=lapply(seq(1:3),function(i){sample.int(100)/10})
    plotCI(data,upperCI,lowerCI,sort=T,stack=T,xlab="Sorted Index",ylab="Value of Random Number:side-by-side",main="Vector Input, Stacked")
    # sine/cosine waves
    data = list();upperCI = list();lowerCI = list()
    data[[1]] = sapply(seq(0,99),sin)
    names(data[[1]]) = seq(0,99)
    upperCI[[1]] = sample(seq(0,1,length.out=1000),100)
    lowerCI[[1]] = sample(seq(0,1,length.out=1000),100)
    data[[2]] = sapply(seq(0,99),cos)
    names(data[[2]]) = seq(0,99)
    upperCI[[2]] = sample(seq(0,1,length.out=1000),100)
    lowerCI[[2]] = sample(seq(0,1,length.out=1000),100)
    plotCI(data,upperCI,lowerCI,sort=F,xlab="x",ylab="sin(x)",main="Multiple Time Series")
    plotCI(data,upperCI,lowerCI,sort=F,stack=T,xlab="x",ylab="sin(x)",main="Multiple Time Series: side-by-side")
  }
  o <- c()
  if (is.null(cols)) cols <- c(brewer.pal(9,"Set3"),brewer.pal(8,"Set3"),brewer.pal(12,"Set3"))
  #try(dev.off(),silent=T)
  if (class(data) == "list") {
    minD <- min(unlist(sapply(data,min)))
    maxD <- max(unlist(sapply(data,max)))
    if (!is.null(lowerCI)) {
      maxCI <- max(sapply(c(upperCI,lowerCI),max))
    } else {
      maxCI <- max(sapply(upperCI,max))
    }
    maxCImin <- minD-maxCI
    maxCImax <- maxD+maxCI
    if (stack) {
      par(mfrow=c(1,length(data)),mar=c(5,2,5,2))
    }
    if (autoX) {
      if (is.null(names(data[[1]]))) {
        xlimits <- sapply(data,length)
        xlimits <- c(1,max(xlimits))
      } else if (length(grep("[A-Za-z]",names(data[[1]])))>0){
        xlimits <- sapply(data,length)
        xlimits <- c(1,max(xlimits))
      } else {
        xlimits <- as.numeric(sapply(data,names))
        xlimits <- c(as.numeric(min(xlimits)),as.numeric(max(xlimits)))
      }
    }
    for (i in 1:length(data)) {
      if (is.null(names(data[[i]]))) {
        # if no x vals are specified, just make into a sequence
        names(data[[i]]) <- seq(1:length(data[[i]]))
      } else if (length(grep("[A-Za-z]",names(data[[i]])))>0){
        names(data[[i]]) <- seq(1:length(data[[i]]))
      }
      if (i == 1) {
        if (sort) {
          o <- order(data[[i]])
          data[[i]] <- data[[i]][order(data[[i]])]
          upperCI[[i]] <- upperCI[[i]][order(data[[i]])]
          if (!is.null(lowerCI)) {
            lowerCI[[i]] <- lowerCI[[i]][order(data[[i]])]
          }
          names(data[[i]]) <- seq(1,length(data[[i]]))
        }
        if (!is.null(NULL)) {
          # Sort data by "time interval"
          data[[i]] = data[[i]][order(as.numeric(names(data[[i]])))]
        }
        # Plot lines
        if (autoX && autoY) {
          plot(names(data[[i]]),data[[i]],col=cols[i],type="l",ylim=c(maxCImin,maxCImax),xlim=xlimits,...)
        } else if (autoY) {
          plot(names(data[[i]]),data[[i]],col=cols[i],type="l",ylim=c(maxCImin,maxCImax),...)
        } else if (autoX) {
          plot(names(data[[i]]),data[[i]],col=cols[i],type="l",xlim=xlimits,...)
        } else {
          plot(names(data[[i]]),data[[i]],col=cols[i],type="l",...)
        }
        # Plot Polygon
        #x <- c(seq(1,length(data[[i]]),1),seq(length(data[[i]]),1,-1))
        x <- c(names(data[[i]]),rev(names(data[[i]])))
        if (!is.null(lowerCI)) {
          y1 <- c(data[[i]]+upperCI[[i]],rev(data[[i]]))
          polygon(x=x,y=y1,col=cols[i],lty=0)
          y2 <- c(data[[i]]-lowerCI[[i]],rev(data[[i]]))
          polygon(x=x,y=y2,col=cols[i],lty=0)
        } else {
          y1 <- c(data[[i]]+upperCI[[i]],rev(data[[i]]))
          polygon(x=x,y=y1,col=cols[i],lty=0)
          y2 <- c(data[[i]]-upperCI[[i]],rev(data[[i]]))
          polygon(x=x,y=y2,col=cols[i],lty=0)
        }
        lines(names(data[[i]]),data[[i]],col=cols[i],type="l",ylim=c(maxCImin,maxCImax),...)
      } else {
        if (sort) { 
          data[[i]] <- data[[i]][o]
          upperCI[[i]] <- upperCI[[i]][o]
          if (!is.null(lowerCI)) {
            lowerCI[[i]] <- lowerCI[[i]][o]
          }
          names(data[[i]]) <- seq(1,length(data[[i]]))
        }
        if (!is.null(NULL)) {
          # Sort data by "time interval"
          data[[i]] = data[[i]][order(as.numeric(names(data[[i]])))]
        }
        if (stack) {
          # Plot lines
          if (autoX && autoY) {
            plot(names(data[[i]]),data[[i]],col=cols[i],type="l",ylim=c(maxCImin,maxCImax),xlim=xlimits,...)
          } else if (autoY) {
            plot(names(data[[i]]),data[[i]],col=cols[i],type="l",ylim=c(maxCImin,maxCImax),...)
          } else if (autoX) {
            plot(names(data[[i]]),data[[i]],col=cols[i],type="l",xlim=xlimits,...)
          } else {
            plot(names(data[[i]]),data[[i]],col=cols[i],type="l",...)
          }
        }
        #x <- c(seq(1,length(data[[i]]),1),seq(length(data[[i]]),1,-1))
        x <- c(names(data[[i]]),rev(names(data[[i]])))
        if (!is.null(lowerCI)) {
          y1 <- c(data[[i]]+upperCI[[i]],rev(data[[i]]))
          polygon(x=x,y=y1,col=cols[i],lty=0)
          y2 <- c(data[[i]]-lowerCI[[i]],rev(data[[i]]))
          polygon(x=x,y=y2,col=cols[i],lty=0)
        } else {
          y1 <- c(data[[i]]+upperCI[[i]],rev(data[[i]]))
          polygon(x=x,y=y1,col=cols[i],lty=0)
          y2 <- c(data[[i]]-upperCI[[i]],rev(data[[i]]))
          polygon(x=x,y=y2,col=cols[i],lty=0)
        }
        lines(names(data[[i]]),data[[i]],col=cols[i],...)
      }
      # Plot legend
      if (legend) {
        legend(position,legend=names(data),fill=cols)
      }
    } 
  } else if (class(data)=="numeric" || class(data)=="integer") {
    minD <- min(data)
    maxD <- max(data)
    if (!is.null(lowerCI)) {
      maxCI <- max(lowerCI)
    } else {
      maxCI <- max(upperCI)
    }
    maxCImin <- minD-maxCI
    maxCImax <- maxD+maxCI
    if (is.null(names(data))) {
        # if no x vals are specified, just make into a sequence
        names(data) <- seq(1:length(data))
    } else if (length(grep("[A-Za-z]",names(data)))>0){
      names(data) <- seq(1:length(data))
    }
    if (sort) {
      o <- order(data)
      data <- data[order(data)]
      upperCI <- upperCI[order(data)]
      if (!is.null(lowerCI)) {
        lowerCI <- lowerCI[order(data)]
      }
      names(data) <- seq(1,length(data))
    }
    if (!is.null(NULL)) {
      # Sort data by "time interval"
      data[[i]] = data[[i]][order(as.numeric(names(data[[i]])))]
    }
    # Plot lines
    if (autoY) {
      plot(names(data),data,col=cols[1],type="l",ylim=c(maxCImin,maxCImax),...)
    } else {
      plot(names(data),data,col=cols[1],type="l",...)
    }
    # Plot Polygon
    #x <- c(seq(1,length(data),1),seq(length(data),1,-1))
    x <- c(names(data),rev(names(data)))
    if (!is.null(lowerCI)) {
      y1 <- c(data+upperCI,rev(data))
      polygon(x=x,y=y1,col=cols[1],lty=0)
      y2 <- c(data-lowerCI,rev(data))
      polygon(x=x,y=y2,col=cols[1],lty=0)
    } else {
      y1 <- c(data+upperCI,rev(data))
      polygon(x=x,y=y1,col=cols[1],lty=0)
      y2 <- c(data-upperCI,rev(data))
      polygon(x=x,y=y2,col=cols[1],lty=0)
    }
    lines(names(data),data,col=cols[1],type="l",ylim=c(maxCImin,maxCImax),...)
  }
  return(o)
}