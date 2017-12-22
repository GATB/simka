


args <- commandArgs(trailingOnly = TRUE)
dir = dirname(sub(".*=", "", commandArgs()[4]))
distanceMatrixFilename = args[1]
#outputFilename = args[2]
pca_axis1 = as.numeric(args[3])
pca_axis2 = as.numeric(args[4])

width = as.numeric(args[5])
height = as.numeric(args[6])
format = args[7]

if(format == "png"){
	png(file=paste0(args[2], ".png"), width=width, height=height, units="in",res=72)
} else{
	pdf(file=paste0(args[2], ".pdf"), width=width, height=height)
}


distanceMatrix = as.matrix(read.table(file=distanceMatrixFilename, sep=";", header=TRUE, row.names=1))
distanceMatrix[lower.tri(distanceMatrix)] <- t(distanceMatrix)[lower.tri(distanceMatrix)] #symmetrize matrix


distance_name = basename(distanceMatrixFilename)
distance_name = unlist(strsplit(distance_name, "[.]"))[1]
distance_name = gsub("mat_", "", distance_name)

use_metadata = F
if(length(args) == 9){
	suppressPackageStartupMessages(library(dendextend))
	
	use_metadata = T
	metadata_table = as.matrix(read.table(file=args[8], sep=";", header=TRUE, row.names=1))
	metadata_variable = args[9]
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

#pdf(file=outputFilename)

par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(x, y, type='n',
xlab=paste0("PC", pca_axis1, " (", round(distData$eig[pca_axis1], 2), "%)"),
ylab=paste0("PC", pca_axis2, " (", round(distData$eig[pca_axis2], 2), "%)")
)

if(use_metadata){
	text(x, y, labels=rownames(distData$data), col=colors_numeric, font=2)
	legend("right", title=metadata_variable, legend=unique(colors), col=unique(colors_numeric), pch=16, inset=c(-0.3, 0))
} else{
	text(x, y, labels=rownames(distData$data), font=2)
}

title(paste0("Simka PCoA\n", distance_name))



