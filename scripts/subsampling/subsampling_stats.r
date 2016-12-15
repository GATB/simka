
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
	
	kendall_correlations = c()
	spearman_correlations = c()
	
	for(filename in distance_matrix_filenames){
		
		distance_matrix = as.matrix(read.table(file=filename, sep=";", header=TRUE, row.names=1))
		distance_matrix = as.dist(distance_matrix)
		
		m = mantel(truth_distance_matrix, distance_matrix, method="spearman", permutations=1)
		correlation = round(m$statistic,3)
		spearman_correlations = c(spearman_correlations, correlation)
		
		
		m = mantel(truth_distance_matrix, distance_matrix, method="kendall", permutations=1)
		correlation = round(m$statistic,3)
		kendall_correlations = c(kendall_correlations, correlation)
	}
	
	return(list(
	spearman_mean=mean(spearman_correlations), spearman_sd=sd(spearman_correlations), spearman_se=sd(spearman_correlations)/sqrt(length(spearman_correlations)),
	kendall_mean=mean(kendall_correlations), kendall_sd=sd(kendall_correlations), kendall_se=sd(kendall_correlations)/sqrt(length(kendall_correlations))
	)
	)
}

#------------------------------------------------------------------

nb_reads = c()
spearman_mean = c()
spearman_sd = c()
spearman_se = c()
kendall_mean = c()
kendall_sd = c()
kendall_se = c()


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
	
	results = compute_mean_correlation(filenames)
	spearman_mean = c(spearman_mean, results$spearman_mean)
	spearman_sd = c(spearman_sd, results$spearman_sd)
	spearman_se = c(spearman_se, results$spearman_se)
	kendall_mean = c(kendall_mean, results$kendall_mean)
	kendall_sd = c(kendall_sd, results$kendall_sd)
	kendall_se = c(kendall_se, results$kendall_se)
	
}
close(inputFile)

#------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create figure
#------------------------------------------------------------------------------------------------------------------------------------------------------------
create_figure_correlation_spearman = function(df){
	
	#print(df)
	#mdf = df
	#mdf <- melt(df, id="ksizes")
	#print(mdf)
	
	p = ggplot(df, aes(x=nb_reads, y=spearman_mean)) +
	geom_point() +
	geom_errorbar(aes(ymin=spearman_mean-spearman_sd, ymax=spearman_mean+spearman_sd), color="green")+
	geom_errorbar(aes(ymin=spearman_mean-spearman_se, ymax=spearman_mean+spearman_se), color="red")+
	ylim(0, 1)

	p = p + theme_bw() +
	theme(
	panel.border = element_rect(colour = "black", fill=NA, size=1),
	text = element_text(size=10)
	)
	p = p + xlab("number of reads")
	p = p + ylab("Spearman correlation")
	
	return(p)
}

create_figure_correlation_kendall = function(df){
	
	print(df)
	#mdf = df
	#mdf <- melt(df, id="ksizes")
	#print(mdf)
	
	p = ggplot(df, aes(x=nb_reads, y=kendall_mean)) +
	geom_point() +
	geom_errorbar(aes(ymin=kendall_mean-kendall_sd, ymax=kendall_mean+kendall_sd), color="green") +
	geom_errorbar(aes(ymin=kendall_mean-kendall_se, ymax=kendall_mean+kendall_se), color="red") +
	ylim(0, 1)

	p = p + theme_bw() +
	theme(
	panel.border = element_rect(colour = "black", fill=NA, size=1),
	text = element_text(size=10)
	)
	p = p + xlab("number of reads")
	p = p + ylab("Kendall correlation")
	
	
	return(p)
}

#------------------------------------------------------------------

df = data.frame(
nb_reads = nb_reads,
spearman_mean = spearman_mean,
spearman_sd = spearman_sd,
spearman_se = spearman_se,
kendall_mean = kendall_mean,
kendall_sd = kendall_sd,
kendall_se = kendall_se
)



p = create_figure_correlation_spearman(df) + ggtitle(paste0("Impact of read filtering - ", distance_name))
ggsave(p, file=paste0("result_figures", "/",  distance_name, "/", "spearman_vs_nbReads.png"))
p = create_figure_correlation_kendall(df) + ggtitle(paste0("Impact of read filtering - ", distance_name))
ggsave(p, file=paste0("result_figures", "/",  distance_name, "/", "kendall_vs_nbReads.png"))
