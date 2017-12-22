# Contributors :
#   Pierre PETERLONGO, pierre.peterlongo@inria.fr [12/06/13]
#   Nicolas MAILLET, nicolas.maillet@inria.fr     [12/06/13]
#   Guillaume Collet, guillaume@gcollet.fr        [27/05/14]
#   Gaetan BENOIT, gaetan.benoit@inria.fr        [08/10/15]
#   Claire LEMAITRE, claire.lemaitre@inria.fr    [06/07/16]
#
# This software is a computer program whose purpose is to find all the
# similar reads between sets of NGS reads. It also provide a similarity
# score between the two samples.
#
# Copyright (C) 2014  INRIA
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.



## Usage : Rscript heatmap.r matrix_asym.csv matrix_sym.csv output_file.pdf title

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
#png(file=args[3],width=800,height=800,res=65)


width = as.numeric(args[4])
height = as.numeric(args[5])
format = args[6]

if(format == "png"){
	png(file=paste0(args[3], ".png"), width=width, height=height, units="in",res=72)
} else{
	pdf(file=paste0(args[3], ".pdf"), width=width, height=height)
}

cr3 = as.matrix(read.table(file=args[1], sep=";", header=TRUE, row.names=1))  # can be symetric matrix
cr3_norm = as.matrix(read.table(file=args[2], sep=";", header=TRUE, row.names=1))  # must be a symetric matrix
cr3[lower.tri(cr3)] <- t(cr3)[lower.tri(cr3)] #symmetrize matrix
cr3_norm[lower.tri(cr3_norm)] <- t(cr3_norm)[lower.tri(cr3_norm)] #symmetrize matrix


distance_name = basename(args[1])
distance_name = unlist(strsplit(distance_name, "[.]"))[1]
distance_name = gsub("mat_", "", distance_name)


use_metadata = F
if(length(args) == 8){
	
	use_metadata = T
	metadata_table = as.matrix(read.table(file=args[7], sep=";", header=TRUE, row.names=1))
	metadata_variable = args[8]
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
	dataset_ids = rownames(cr3_norm)
	for(i in 1:dim(cr3_norm)[1]){
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





n=100 # number of steps between 2 colors

## Transforming 0-1 distances in 0-100 similarity measure
if(grepl("chord",args[1]) || grepl("hellinger",args[1])){
	cr3 = (sqrt(2) - cr3) * 100
} else {
	cr3 = (1 - cr3) * 100
}


## Computing mini-maxi for colour palette

mini=min(cr3[])
maxi=max(cr3[row(cr3)!=col(cr3)]) # ignoring the diagonal
trueMax=max(cr3[]) # typically the value in the diagonal = 100
q25=quantile(cr3[row(cr3)!=col(cr3)],0.25,1)
q50=quantile(cr3[row(cr3)!=col(cr3)],0.5,1)
q75=quantile(cr3[row(cr3)!=col(cr3)],0.75,1)

## We use the quantiles to ignore some outlier values in the matrix (values<mini will have colour of mini and values>maxi will have a colour between brown and grey23)
mini=max(q25-1.5*(q75-q25),0)
maxi=min(q75+1.5*(q75-q25),trueMax)

palette=colorRampPalette(c("green", "yellow", "red", "brown", "grey23"))(n = 5*n-1)

## Checking if maxi = trueMax
trueMax.needed=ifelse(maxi<trueMax,"T","F")

if(trueMax.needed){
  breaks=c(seq(mini,maxi,length=4*n),seq(maxi+1e-5,trueMax,length=n))
  # breaks are equally distributed in the range mini-maxi (intervals can be different in the range maxi-trueMax, containing very few points)
} else {
  breaks=c(seq(mini,maxi,length=5*n))
}



# Dendrogram is obtained with the symetric matrix
 distance    = dist(cr3_norm)
 cluster     = hclust(distance, method="average")
dendrogram  = as.dendrogram(cluster)

# Heatmap 
par(fig=c(0.2,1,0,0.8),mar=rep(1,4))

if(use_metadata){
	
	colors_numeric=as.character(colors_numeric)
 heatmap.2(cr3,
 trace = "none",
 dendrogram = "row",
 key = FALSE,
 Rowv=dendrogram,
 Colv = rev(dendrogram),
 col=palette,
 breaks = breaks,
 margins=c(10,10),
 main=paste0("Simka heatmap\n", distance_name), sub="", cexRow = 0.8, cexCol = 0.8, RowSideColors=colors_numeric, ColSideColors=colors_numeric)
	
} else {
	
 heatmap.2(cr3,
 trace = "none",
 dendrogram = "row",
 key = FALSE,
 Rowv=dendrogram,
 Colv = rev(dendrogram),
 col=palette,
 breaks = breaks,
 margins=c(10,10),
 main=paste0("Simka heatmap\n", distance_name), sub="", cexRow = 0.8, cexCol = 0.8)
}



if(use_metadata){
	par(lend = 1)           # square line ends for the color legend
	legend("topright", title=metadata_variable, legend=unique(colors), col=unique(colors_numeric), pch=16, lty= 1,             # line style
	lwd = 10
	)
	
}



# Adding the colour scale
par(fig=c(0.05,0.4,0.8,1), mar=rep(2,4), new=TRUE)

if(trueMax.needed){

  diff=maxi-mini
  breaksToMaxi=breaks[1:(4*n)] # using only breaks from mini to maxi
  black.width=max(diff/9)
  black.space=max(diff/9)
  
  plot(c(mini,maxi+black.width+black.space),c(0,2),type="n",yaxt="n",ylab="",xlab="",xaxt="n",xaxs = "i", yaxs = "i")
  rect(breaksToMaxi[-length(breaksToMaxi)],0,breaksToMaxi[-1],2,col=palette,border=NA)

  
  ti=pretty(breaksToMaxi)
  ti=ti[ti<maxi]
  axis(1,at=c(ti,maxi+black.space+black.width/2),label=c(ti,trueMax))
  
  # Here plotting the TrueMax colour with a white space
  rect(maxi+black.space,0,maxi+black.space+black.width,2,col=palette[5*n-1],border=NA)
  rect(maxi,-0.1,maxi+black.space,2.1,col="white",border=NA)

} else{
 plot(range(breaks),c(0,2),type="n",yaxt="n",ylab="",xlab="",xaxs = "i", yaxs = "i")
 rect(breaks[-length(breaks)],0,breaks[-1],2,col=palette,border=NA)
}






d=dev.off()


