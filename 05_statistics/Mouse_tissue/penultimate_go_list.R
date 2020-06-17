
######################################################
######################################################
## This script is used to process penultimate exons ##
######################################################
######################################################

library(reshape2)
library(ggplot2)
library(ggforce)
library(textshape)
library(dplyr)
library(sva)
library(weights)

##################################
##### prepare PSI data frame #####
##################################

## load tissue type
meta <- read.table("/scratch_space/qgao/Penultimate/2.data/mouse_SRR_Acc_List.txt", stringsAsFactor=F, row.names = 1)
colnames(meta) <- c("Tissue", "Library")

## load data
penul <- read.table("penultimate_psi.txt", header = T, stringsAsFactor=F)

## PSI
psi <- penul[, c(1, 4:112)]

###########################################
##### generate RNA-seq support figure #####
###########################################

## melt
penul_melt <- melt(psi)

## add tissue
penul_melt$Tissue <- meta[match(penul_melt$variable, rownames(meta)), 1]

## split by tissue
penul_split <- split(penul_melt, penul_melt$Tissue)

## calculate average psi for each tissue
penul_avg_list = list()
for(i in 1:length(penul_split))
{
	tmp = split(penul_split[[i]], penul_split[[i]][,1])
	penul_avg_list[[i]] = lapply(tmp, function(x) mean(x$value, na.rm=T))
}
penul_avg <- do.call(cbind, penul_avg_list)
colnames(penul_avg) <- sort(unique(meta$Tissue))
penul_avg[penul_avg=="NaN"] <- 1000

## load gene id
intron <- read.table("/scratch_space/qgao/Penultimate/1.reference/mouse/F_MISO/Eligibile/mm10.SE.intron.bed", stringsAsFactor=F)
intron$Gene <- matrix(unlist(strsplit(intron[,4], "__")), byrow=T, ncol=3)[,2]
event_gene <- intron[seq(1, 6466, 2), c(5,7)]

## add gene
penul_tissue <- as.data.frame(penul_avg)
penul_tissue$Gene <- event_gene[match(rownames(penul_tissue), event_gene[,1]),2]

## get number of events and genes for supp table
ncount <- vector()
## output
for(i in 1:22)
{
	tmp <- penul_tissue[penul_tissue[,i]<90, 23, drop=F]
	write.table(tmp, paste0(colnames(penul_tissue)[i], ".gene.tsv"), row.names=T, col.names=F, sep="\t", quote=F)
	ncount <- rbind(ncount, c(colnames(penul_tissue)[i], nrow(tmp), length(unique(tmp[,1]))))
}

## get all
suppall <- penul_tissue[apply(penul_tissue[,1:22]<90, 1, all), 23, drop=F]
write.table(suppall, "alltissue.gene.tsv", row.names=T, col.names=F, sep="\t", quote=F)
ncount <- rbind(ncount, c("All", nrow(suppall), length(unique(suppall[,1]))))

## get any
suppany <- penul_tissue[apply(penul_tissue[,1:22]<90, 1, any), 23, drop=F]
write.table(suppany, "anytissue.gene.tsv", row.names=T, col.names=F, sep="\t", quote=F)
ncount <- rbind(ncount, c("Any", nrow(suppany), length(unique(suppany[,1]))))

write.table(ncount, "Numbers_for_GO.txt", quote=F, row.names=F, col.names=F, sep="\t")

