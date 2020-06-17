
###################################################
###################################################
## This script is used to process internal exons ##
###################################################
###################################################

library(reshape2)
library(ggplot2)

##################################
##### prepare PSI data frame #####
##################################

## load tissue type
meta <- read.table("/scratch_space/qgao/Penultimate/2.data/paired_human_SRR_Acc_List.txt", stringsAsFactor=F, row.names = 1)
colnames(meta) <- c("Tissue")

## load data
internal <- read.table("pen_3n.txt", header = T, stringsAsFactor=F)

## PSI
psi <- internal[, c(1, 3:18)]

###########################################
##### generate RNA-seq support figure #####
###########################################

## melt
internal_melt <- melt(psi)

## add tissue
internal_melt$Tissue <- meta[match(internal_melt$variable, rownames(meta)), 1]

## split by tissue
internal_split <- split(internal_melt, internal_melt$Tissue)

## calculate average psi for each tissue
internal_avg_list = list()
for(i in 1:length(internal_split))
{
	tmp = split(internal_split[[i]], internal_split[[i]][,1])
	internal_avg_list[[i]] = lapply(tmp, function(x) mean(x$value, na.rm=T))
}
internal_avg <- do.call(cbind, internal_avg_list)
colnames(internal_avg) <- sort(unique(meta$Tissue))
internal_avg[internal_avg=="NaN"] <- 1000

## calculate number of alternatively spliced exons
internal_as <- data.frame(Number=apply(internal_avg<90, 1, sum))
internal_as$Annotation <- internal[match(rownames(internal_as), internal$Event),2]
internal_as$Type <- internal[match(rownames(internal_as), internal$Event),3]
internal_as$Support <- ifelse(internal_as$Number > 0, 1, 0)
## statistics
length(internal_as$Support)
#total: 9906
sum(internal_as$Support)
#supported: 6377

## output
write.table(internal_as, "human_3N_support_info.txt", row.names=T, col.names=NA, quote=F, sep="\t")

## calculate number of alternatively spliced exons
## cutoff 85
internal_as <- data.frame(Number=apply(internal_avg<85, 1, sum))
internal_as$Annotation <- internal[match(rownames(internal_as), internal$Event),2]
internal_as$Type <- internal[match(rownames(internal_as), internal$Event),3]
internal_as$Support <- ifelse(internal_as$Number > 0, 1, 0)
sum(internal_as$Support)
#supported: 4771

## calculate number of alternatively spliced exons
## cutoff 80
internal_as <- data.frame(Number=apply(internal_avg<80, 1, sum))
internal_as$Annotation <- internal[match(rownames(internal_as), internal$Event),2]
internal_as$Type <- internal[match(rownames(internal_as), internal$Event),3]
internal_as$Support <- ifelse(internal_as$Number > 0, 1, 0)
sum(internal_as$Support)
#supported: 4110


