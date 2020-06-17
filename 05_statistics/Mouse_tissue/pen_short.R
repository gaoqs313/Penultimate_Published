
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
meta <- read.table("/scratch_space/qgao/Penultimate/2.data/mouse_SRR_Acc_List.txt", stringsAsFactor=F, row.names = 1)
colnames(meta) <- c("Tissue", "Library")

## load data
internal <- read.table("pen_short.txt", header = T, stringsAsFactor=F)

## PSI
psi <- internal[, c(1, 3:111)]

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
#total: 9791
sum(internal_as$Support)
#supported: 6817

## output
write.table(internal_as, "mouse_short_support_info.txt", row.names=T, col.names=NA, quote=F, sep="\t")

## calculate number of alternatively spliced exons
## cutoff 85
internal_as <- data.frame(Number=apply(internal_avg<85, 1, sum))
internal_as$Annotation <- internal[match(rownames(internal_as), internal$Event),2]
internal_as$Type <- internal[match(rownames(internal_as), internal$Event),3]
internal_as$Support <- ifelse(internal_as$Number > 0, 1, 0)
sum(internal_as$Support)
#supported: 5712

## calculate number of alternatively spliced exons
## cutoff 80
internal_as <- data.frame(Number=apply(internal_avg<80, 1, sum))
internal_as$Annotation <- internal[match(rownames(internal_as), internal$Event),2]
internal_as$Type <- internal[match(rownames(internal_as), internal$Event),3]
internal_as$Support <- ifelse(internal_as$Number > 0, 1, 0)
sum(internal_as$Support)
#supported: 4955

###############################################################
## calculate mean and maximum PSI distribution
penul_avg <- internal_avg
penul <- internal
penul_as <- data.frame(Number=apply(penul_avg<90, 1, sum))
penul_as$Annotation <- penul[match(rownames(penul_as), penul$Event),2]
penul_as$Type <- penul[match(rownames(penul_as), penul$Event),3]
penul_as$Support <- ifelse(penul_as$Number > 0, 1, 0)

penul_tissue <- as.data.frame(penul_avg)
penul_tissue <- sapply(penul_tissue, as.numeric)
rownames(penul_tissue) <- rownames(penul_avg)

## get RNA-seq supported events
penul_tissue_sup <- penul_tissue[rownames(penul_tissue) %in% rownames(penul_as[penul_as$Support==1,]), ]

## calculate mean and maximum PSI
penul_tissue_val <- data.frame(PSI_MIN=apply(penul_tissue_sup, 1, function(x) min(x[x<1000])), PSI_MEDIAN=apply(penul_tissue_sup, 1, function(x) median(x[x<1000])), PSI_MEAN=apply(penul_tissue_sup, 1, function(x) mean(x[x<1000])), PSI_MAX=apply(penul_tissue_sup, 1, function(x) max(x[x<1000])))

## make a copy for combined set
penul_tissue_val2 <- penul_tissue_val
penul_tissue_val2$Annotation <- "All"

## add annotation
penul_tissue_val$Annotation <- penul[match(rownames(penul_tissue_val), penul$Event),2]
penul_tissue_val <- rbind(penul_tissue_val, penul_tissue_val2)

## melt data
df_psi <- melt(penul_tissue_val)

## plot
pdf("Aggregrated_PSI_Short_distribution.pdf", width=8)
p <- ggplot(df_psi, aes(x=Annotation, y=value, fill=Annotation)) +
  geom_violin() +
  geom_boxplot(width=0.1, outlier.color=NA)+
  theme_bw() +
  facet_wrap(.~variable) +
  scale_fill_manual(name="", values=c("#B3B3B3", "#A6D854", "#E78AC3"), labels=c("All supported (6817)", "Annotated (419, 6.1%)", "Novel (6398, 93.9%)")) +
  xlab("") +
  ylab("Percent Splice In") +
  theme(legend.key = element_blank(), strip.background = element_rect(fill="white"))+
  theme(axis.text.x = element_text(size = 10, family = "Helvetica")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=18, family = "Helvetica")) +
  theme(legend.position = "top")
print(p)
dev.off()


