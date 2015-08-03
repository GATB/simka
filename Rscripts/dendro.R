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


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
png(file=args[2],width=800,height=800,res=100)
cr3 = as.matrix(read.table(file=args[1], sep=";", header=TRUE, row.names=1))
inv_cr3 = matrix(100, ncol=dim(cr3)[1], nrow=dim(cr3)[1]) - cr3
Commet_distance = as.dist(inv_cr3)
dendo_cr3 = hclust(Commet_distance)
plot(dendo_cr3, main="Commet normalized analysis", sub = NA, xlab = paste("Complete clusterization of ",args[1]))
