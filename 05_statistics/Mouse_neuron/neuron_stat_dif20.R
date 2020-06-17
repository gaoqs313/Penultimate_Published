
######################################################
######################################################
## This script is used to process penultimate exons ##
######################################################
######################################################

library(reshape2)
library(ggplot2)
library(dplyr)
library(weights)
library(Rtsne)
library(RColorBrewer)

##################################
##### prepare PSI data frame #####
##################################

## load tissue type
meta <- read.table("/scratch_space/qgao/Penultimate/2.data/neuron_SRR_Acc_List.txt", stringsAsFactor=F)
colnames(meta) <- c("Stage", "Tissue")

## load data
psi <- read.table("/scratch_space/qgao/Penultimate/11.revision/E_Mouse_GO/Neuron/neuron_psi.txt", header = T, stringsAsFactor=F)

###########################################
##### generate RNA-seq support figure #####
###########################################

## melt
penul_melt <- melt(psi)

## add tissue
penul_melt$Tissue <- meta[match(penul_melt$variable, meta[,1]), 2]

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

## calculate number of alternatively spliced exons
penul_as <- penul_avg[apply(penul_avg < 90, 1, any), ]
penul_as <- apply(penul_as, 2, unlist)
#775

## get differential events
penul_as <- penul_as[apply(penul_as<1000, 1, all), ]
penul_dif <- penul_as[apply(penul_as, 1, function(x) max(x[x<1000])-min(x[x<1000])) > 20, ] 
#276

## load gene id
intron <- read.table("/scratch_space/qgao/Penultimate/1.reference/mouse/F_MISO/Eligibile/mm10.SE.intron.bed", stringsAsFactor=F)
intron$Gene <- matrix(unlist(strsplit(intron[,4], "__")), byrow=T, ncol=3)[,2]
event_gene <- intron[seq(1, 6466, 2), c(5,7)]

## add gene
penul_tissue <- as.data.frame(penul_dif)
penul_tissue$Gene <- event_gene[match(rownames(penul_tissue), event_gene[,1]),2]

## output
output <- penul_tissue[, c(2,1,3,4,8,5:7,9)]
write.table(output, "neuron_diff.tsv", row.names=T, col.names=NA, quote=F, sep="\t")

## background
penul_con <- penul_avg[apply(penul_avg<1000, 1, all), ]
penul_con <- penul_con[apply(penul_con>90, 1, all) | apply(penul_con<30, 1, all),]
penul_con <- apply(penul_con, 2, unlist)
#838 + 57 = 895

output2 <- penul_con[, c(2,1,3,4,8,5:7)]
write.table(output2, "neuron_same.tsv", row.names=T, col.names=NA, quote=F, sep="\t")
