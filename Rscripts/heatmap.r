# Contributors :
#   Pierre PETERLONGO, pierre.peterlongo@inria.fr [12/06/13]
#   Nicolas MAILLET, nicolas.maillet@inria.fr     [12/06/13]
#   Guillaume Collet, guillaume@gcollet.fr        [27/05/14]
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

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
png(file=args[3],width=800,height=800,res=65)
n=100 # number of steps between 2 colors
cr3 = as.matrix(read.table(file=args[1], sep=";", header=TRUE, row.names=1))

mini=min(cr3[])
maxi=max(cr3[row(cr3)!=col(cr3)])
trueMax=max(cr3[])
q25=quantile(cr3[row(cr3)!=col(cr3)],0.25,1)
q50=quantile(cr3[row(cr3)!=col(cr3)],0.5,1)
q75=quantile(cr3[row(cr3)!=col(cr3)],0.75,1)

mini=max(q25-1.5*(q75-q25),0)
maxi=min(q75+1.5*(q75-q25),trueMax)
diff=maxi-mini

palette=colorRampPalette(c("green", "yellow", "red", "brown", "grey23"))(n = 5*n-1)

 breaks=c(seq(mini,mini+diff/4,length=n), # for green
                seq(mini+diff/4,mini+diff/2,length=n), # for yellow
                seq(mini+diff/2,mini+3*diff/4,length=n), # for red
                seq(mini+3*diff/4,maxi,length=n), # for brown
                seq(maxi,trueMax,length=n)) # for black

 par(fig=c(0.2,1,0,0.8))

 
 cr3_norm = as.matrix(read.table(file=args[2], sep=";", header=TRUE, row.names=1))
 inv_cr3 = matrix(trueMax, ncol=dim(cr3_norm)[1], nrow=dim(cr3_norm)[1]) - cr3_norm
 distance    = dist(inv_cr3)
 cluster     = hclust(distance)
 dendrogram  = as.dendrogram(cluster)
 
 
 heatmap.2(cr3,
 trace = "none",
 dendrogram = "none",
 key = FALSE,
  Rowv=dendrogram,
  Colv = rev(dendrogram),
 col=palette,
 breaks = breaks,
 margins=c(10,10),
 main = args[4])

par(fig=c(0.05,0.4,0.8,1), new=TRUE)

  trueMin=min(cr3[])

 breaksToMaxi=breaks[1:(4*n)] # prend que les breaks <=maxi
 black.width=max(diff/9)
 black.space=max(diff/9)


 par(xpd=T,cex=1.8,mar=c(1,1,1,1))
 plot(c(trueMin,maxi+black.width+black.space),c(0,2),type="n",yaxt="n",ylab="",xlab="",xaxt="n",xaxs = "i", yaxs = "i")
 rect(breaksToMaxi[-length(breaksToMaxi)],0,breaksToMaxi[-1],2,col=palette,border=NA)
 rect(maxi+black.space,0,maxi+black.space+black.width,2,col=palette[5*n-1],border=NA)

 ti=pretty(breaksToMaxi)
 ti=ti[ti<maxi]
 trueMax
 axis(1,at=c(ti,maxi+black.space+black.width/2),label=c(ti,trueMax))

 # pour faire un break
 rect(maxi,-0.1,maxi+black.space,2.1,col="white",border=NA)





