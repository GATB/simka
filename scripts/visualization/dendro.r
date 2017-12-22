#Author: Gaetan Benoit
#Contact: gaetan.benoit@inria.fr


args <- commandArgs(trailingOnly = TRUE)
distanceMatrixFilename = args[1]
distance_name = basename(distanceMatrixFilename)
distance_name = unlist(strsplit(distance_name, "[.]"))[1]
distance_name = gsub("mat_", "", distance_name)


distanceMatrix = as.matrix(read.table(file=distanceMatrixFilename, sep=";", header=TRUE, row.names=1))
distanceMatrix[lower.tri(distanceMatrix)] <- t(distanceMatrix)[lower.tri(distanceMatrix)] #symmetrize matrix


width = as.numeric(args[3])
height = as.numeric(args[4])
format = args[5]

if(format == "png"){
	png(file=paste0(args[2], ".png"), width=width, height=height, units="in",res=72)
} else{
	pdf(file=paste0(args[2], ".pdf"), width=width, height=height)
}


use_metadata = F
if(length(args) == 7){
	suppressPackageStartupMessages(library(dendextend))
	
	use_metadata = T
	metadata_table = as.matrix(read.table(file=args[6], sep=";", header=TRUE, row.names=1))
	metadata_variable = args[7]
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
		print(paste0(dataset_id, " ", variables[[i]]))
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




distanceMatrix = distanceMatrix*100
#inv_cr3 = matrix(100, ncol=dim(cr3)[1], nrow=dim(cr3)[1]) - cr3
Commet_distance = as.dist(distanceMatrix)
hc = hclust(Commet_distance, method="average")
dendo_cr3 = as.dendrogram(hc)

if(use_metadata){
	
	colors_numeric_hc = colors_numeric[hc$order]
	dendo_cr3 %>% set("labels_col", colors_numeric_hc) %>% set("branches_k_color", colors_numeric_hc) %>% # change color
	plot(main=paste0("Simka hierarchical clustering\n", distance_name), cex = 0.3, xlab="", sub="")
	legend("topright", title=metadata_variable, legend=unique(colors), col=unique(colors_numeric), pch=16)

} else{
	plot(dendo_cr3, main=paste0("Simka hierarchical clustering\n", distance_name), cex = 0.3, xlab="", sub="")

}




