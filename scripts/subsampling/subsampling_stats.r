
library(ggplot2)
library(grid)
library(gridExtra)
library("reshape2")
require(gtable)
library(vegan)


args <- commandArgs(trailingOnly = TRUE)

distance_name = args[1]
input_filename = args[2]
output_dir_temp = args[3]

setwd(output_dir_temp)


#------------------------------------------------------------------
truth_distance_matrix_filename = paste0(output_dir_temp, "/truth_results/", "mat_", distance_name, ".csv.gz")
print(truth_distance_matrix_filename)
truth_distance_matrix = as.matrix(read.table(file=truth_distance_matrix_filename, sep=";", header=TRUE, row.names=1))
truth_distance_matrix = as.dist(truth_distance_matrix)
#print(truth_distance_matrix)
#distanceMatrix[] <- as.numeric(distanceMatrix)
#------------------------------------------------------------------

compute_mean_correlation = function(distance_matrix_filenames){
	
	correlations = c()
	
	for(filename in distance_matrix_filenames){
		
		distance_matrix = as.matrix(read.table(file=filename, sep=";", header=TRUE, row.names=1))
		distance_matrix = as.dist(distance_matrix)
		
		m = mantel(truth_distance_matrix, distance_matrix, method="spearman", permutations=1)
		correlation = round(m$statistic,3)
		correlations = c(correlations, correlation)
	}
	
	return(mean(correlations))
}

#------------------------------------------------------------------

nb_reads = c()
mean_correlations = c()


inputFile = file(input_filename,open="r")
for(line in readLines(inputFile)){
	fields = strsplit(line, " ")
	fields = fields[[1]]
	
	
	nbReads = as.numeric(fields[1])
	nb_reads = c(nb_reads, nbReads)
	
	filenames = c()
	for(i in 2:length(fields)){
		filenames = c(filenames, fields[i])
	}
	
	mean_correlation = compute_mean_correlation(filenames)
	mean_correlations = c(mean_correlations, mean_correlation)
	
}
close(inputFile)

#------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create figure
#------------------------------------------------------------------------------------------------------------------------------------------------------------
create_figure_correlation = function(df){
	
	#mdf = df
	#mdf <- melt(df, id="ksizes")
	#print(mdf)
	
	p = ggplot(df, aes(x=nb_reads, y=mean_correlations)) +
	geom_point() +
	geom_line()
	
	p = p + theme_bw() +
	theme(
	panel.border = element_rect(colour = "black", fill=NA, size=1),
	text = element_text(size=10)
	)
	p = p + ggtitle("Functional correlation")
	#p = p + scale_colour_manual(values=COLOR_PALETTE)
	
	return(p)
	
}

#------------------------------------------------------------------

df = data.frame(
nb_reads = nb_reads,
mean_correlations = mean_correlations
)


p = create_figure_correlation(df)
ggsave(p, file="result_figures/correlation_vs_nbReads.png")

