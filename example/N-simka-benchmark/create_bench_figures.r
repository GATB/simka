
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library("reshape2")

dir = dirname(sub(".*=", "", commandArgs()[4]))






###################################################################################################
# Time
###################################################################################################
create_time_figure = function(data){
	
	#mdf <- melt(data, id.vars="nb", value.name="value", variable.name="Year")
	mdf <- melt(data, id = "nb_datasets")
	mdf = mdf[!is.na(mdf$value),]
	print(mdf)
	
	p = ggplot(mdf, aes(x=nb_datasets, y=value, group=variable, color=variable)) +
	geom_point() +
	geom_line() +

#scale_x_continuous(breaks=x_scale) +
#scale_y_continuous(expand = c(0, 0)) +
	theme_bw() +
	theme(axis.line = element_line(colour = "black"),
	panel.border = element_blank(),
	panel.background = element_blank(),
	legend.title=element_blank(),
	legend.key = element_blank()
	#legend.margin=unit(0,"cm"),
	#legend.margin=unit(-2, "cm"),
	#plot.margin = unit(x = c(0, 0, 0, 0), units = "cm")
	) +
	xlab("Number of datasets")+
	ylab("CPU time (s)")
	
	return(p)
}

data = data.frame(read.table(paste(dir, "bench_total_time.csv", sep="/"), sep=";", header=TRUE))
print(data)
p = create_time_figure(data) + ggtitle("Total simka time")
ggsave(p, file=paste(dir, "plot_total_time.png", sep="/"), width=6, height=6)

data = data.frame(read.table(paste(dir, "bench_pass_time.csv", sep="/"), sep=";", header=TRUE))
print(data)
p = create_time_figure(data) + ggtitle("Total simka time (non cumulated)")
ggsave(p, file=paste(dir, "plot_pass_time.png", sep="/"), width=6, height=6)

