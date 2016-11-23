#library(ape)

library(ggplot2)
library(grid)
library(gridExtra)
library("reshape2")
require(gtable)
library(vegan)
library(scales)


args <- commandArgs(trailingOnly = TRUE)
dir = dirname(sub(".*=", "", commandArgs()[4]))
distanceMatrixFilename = args[1]


distanceMatrix = as.matrix(read.table(file=distanceMatrixFilename, sep=";", header=TRUE, row.names=1))
distanceMatrix[] <- as.numeric(distanceMatrix)


## ----pca-function--------------------------------------------------------
## Function for plotting nice PCA
create_pcoa_figure <- function(axes = c(1, 2), data, eig) {
	#group = "Enterotype"
	axisName = "Dim"
	sep = ""
	
	#print(data)
	p <- ggplot(data, aes(x = Dim1, y = Dim2, label=rownames(data))) + theme_bw()
	
	#p <- p + geom_point(aes(fill=abundance_vector), colour="black",pch=21, size=2) #, alpha=0.7)
	#p <- p + geom_point(size=1) #, alpha=0.7)
	p <- p + geom_text() #, alpha=0.7)
	
	p <- p + xlab(paste0("PC", axes[1], " (", round(eig[axes[1]], 2), "%)"))
	p <- p + ylab(paste0("PC", axes[2], " (", round(eig[axes[2]], 2), "%)")) +
	
	theme(
	axis.line=element_blank(),axis.text.x=element_blank(),
	axis.text.y=element_blank(),axis.ticks=element_blank(),
	panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(),
	#panel.background = element_rect(fill="#DDDDDD"),
	panel.border = element_rect(colour = "black", fill=NA, size=1)
	#legend.position="none"
	) +
	ggtitle("Simka MDS")
	
	return(p)
}

format_distance <- function(distanceMatrix) {
	distances <- as.dist(distanceMatrix)
	#distances <- as.dist(read.csv2(listDistances[name], row.names = 1))
	numberAxis <- 2
	pcoa <- cmdscale(distances, k = numberAxis, eig = TRUE)
	#print(pcoa)
	plotData <- data.frame(pcoa$points)
	relEig <- (100 * pcoa$eig / sum(abs(pcoa$eig)) )[1:numberAxis]
	colnames(plotData) <- paste0("Dim", seq(1, numberAxis))
	plotData$ID <- rownames(plotData)
	#plotData$Enterotype <- bodysites #ent #gsub(".*_ET", "ET", plotData$ID)
	#print(plotData)
	return(list(data = plotData, eig = relEig))
}

#print(distanceMatrix)

distData <- format_distance(distanceMatrix)
p <- create_pcoa_figure(c(1, 2), distData$data, distData$eig)

ggsave(p, file=args[2])



