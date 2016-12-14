


args <- commandArgs(trailingOnly = TRUE)
dir = dirname(sub(".*=", "", commandArgs()[4]))
distanceMatrixFilename = args[1]
outputFilename = args[2]
pca_axis1 = as.numeric(args[3])
pca_axis2 = as.numeric(args[4])

distanceMatrix = as.matrix(read.table(file=distanceMatrixFilename, sep=";", header=TRUE, row.names=1))

use_metadata = F
if(length(args) == 6){
	suppressPackageStartupMessages(library(dendextend))
	
	use_metadata = T
	metadata_table = as.matrix(read.table(file=args[5], sep=";", header=TRUE, row.names=1))
	metadata_variable = args[6]
	#print(metadata_table)
	variables = metadata_table[,metadata_variable]
	#print(variables)
	
	meatadata_index = list()
	dataset_ids = rownames(metadata_table)
	for(i in 1:length(dataset_ids)){
		dataset_id = dataset_ids[i]
		#print(dataset_id)
		#print(variables[[i]])
		meatadata_index[[dataset_id]] = variables[[i]]
		#print(meatadata_index[[dataset_id]])
	}
	
	colors = c()
	dataset_ids = rownames(distanceMatrix)
	for(i in 1:dim(distanceMatrix)[1]){
		dataset_id = dataset_ids[i]
		colors = c(colors, meatadata_index[[dataset_id]])
	}
	colors_numeric_temp = c()
	colors_numeric = as.numeric(as.factor(colors))
	for(i in 1:length(colors_numeric)){
		colors_numeric_temp = c(colors_numeric_temp, colors_numeric[i]+1)
	}
	colors_numeric = colors_numeric_temp
	#print(colors)
}



#print(metadata_table)


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
distData <- format_distance(distanceMatrix, max(pca_axis1, pca_axis2))
x = distData$data[,paste0("Dim", pca_axis1)]
y = distData$data[,paste0("Dim", pca_axis2)]

pdf(file=outputFilename)

plot(x, y, type='n',
xlab=paste0("PC", pca_axis1, " (", round(distData$eig[pca_axis1], 2), "%)"),
ylab=paste0("PC", pca_axis2, " (", round(distData$eig[pca_axis2], 2), "%)")
)

if(use_metadata){
	text(x, y, labels=rownames(distData$data), col=colors_numeric, font=2)
	legend("right", title=metadata_variable, legend=unique(colors), col=unique(colors_numeric), pch=16)
} else{
	text(x, y, labels=rownames(distData$data), font=2)
}

title("Simka MDS/PCoA")



