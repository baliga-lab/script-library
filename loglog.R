# to plot log-log plot

require(lattice)
require(latticeExtra)

xyplot(y~x, as.data.frame(data), scales=list(y=list(log=10),x=list(log=10)), yscale.components=yscale.components.log10ticks,xscale.components=xscale.components.log10ticks)