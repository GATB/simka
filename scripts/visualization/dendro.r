#Author: Gaetan Benoit
#Contact: gaetan.benoit@inria.fr


args <- commandArgs(trailingOnly = TRUE)
distanceMatrixFilename = args[1]


distanceMatrix = as.matrix(read.table(file=distanceMatrixFilename, sep=";", header=TRUE, row.names=1))

pdf(file=args[2])

print(length(args))

use_metadata = F
if(length(args) == 4){
	suppressPackageStartupMessages(library(dendextend))
	
	use_metadata = T
	metadata_table = as.matrix(read.table(file=args[3], sep=";", header=TRUE, row.names=1))
	metadata_variable = args[4]
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




distanceMatrix = distanceMatrix*100
#inv_cr3 = matrix(100, ncol=dim(cr3)[1], nrow=dim(cr3)[1]) - cr3
Commet_distance = as.dist(distanceMatrix)
hc = hclust(Commet_distance, method="ward.D2")
dendo_cr3 = as.dendrogram(hc)

if(use_metadata){
	
	colors_numeric = colors_numeric[hc$order]
	dendo_cr3 %>% set("labels_col", colors_numeric) %>% set("branches_k_color", colors_numeric) %>% # change color
	plot(main="Simka hierarchical clustering", cex = 0.3, xlab="", sub="")
	legend("topright", title=metadata_variable, legend=unique(colors), col=unique(colors_numeric), pch=16)

} else{
	plot(dendo_cr3, main="Simka hierarchical clustering", cex = 0.3, xlab="", sub="")

}




