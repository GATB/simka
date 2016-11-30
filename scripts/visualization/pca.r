

PCA_AXIS1 = 1
PCA_AXIS2 = 2


#------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
dir = dirname(sub(".*=", "", commandArgs()[4]))
distanceMatrixFilename = args[1]


distanceMatrix = as.matrix(read.table(file=distanceMatrixFilename, sep=";", header=TRUE, row.names=1))
distanceMatrix[] <- as.numeric(distanceMatrix)


format_distance <- function(distanceMatrix, numberAxis) {
	distances <- as.dist(distanceMatrix)
	#distances <- as.dist(read.csv2(listDistances[name], row.names = 1))
	#numberAxis <- 2
	pcoa <- cmdscale(distances, k = numberAxis, eig = TRUE)
	#print(pcoa)
	plotData <- data.frame(pcoa$points)
	relEig <- (100 * pcoa$eig / sum(abs(pcoa$eig)) )[1:numberAxis]
	colnames(plotData) <- paste0("Dim", seq(1, numberAxis))
	#plotData$ID <- rownames(plotData)
	#plotData$Enterotype <- bodysites #ent #gsub(".*_ET", "ET", plotData$ID)
	#print(plotData)
	return(list(data = plotData, eig = relEig))
}

#print(distanceMatrix)
distData <- format_distance(distanceMatrix, max(PCA_AXIS1, PCA_AXIS2))
x = distData$data[,paste0("Dim", PCA_AXIS1)]
y = distData$data[,paste0("Dim", PCA_AXIS2)]

pdf(file=args[2])
plot(x, y, type='n',
xlab=paste0("PC", PCA_AXIS1, " (", round(distData$eig[PCA_AXIS1], 2), "%)"),
ylab=paste0("PC", PCA_AXIS2, " (", round(distData$eig[PCA_AXIS2], 2), "%)"),
)
text(x, y, labels=rownames(distData$data))
title("Simka MDS")



