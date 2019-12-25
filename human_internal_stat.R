
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
meta <- read.table("paired_human_SRR_Acc_List.txt", stringsAsFactor=F, row.names = 1)
colnames(meta) <- c("Tissue")

## load data
internal <- read.table("human_internal_psi.txt", header = T, stringsAsFactor=F)

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
#13805 total
sum(internal_as$Support)
#8219 supported
table(internal_as[,c(2,4)])
#              Support
#Annotation     0    1
# Annotated   395  790
#     Novel  7875 7508

## compare penultimate with internal
bar_data <- as.data.frame(matrix(c(
  paste("Internal","(1185)", sep="\n"), "Annotated", "No", "395", "33.3",
  paste("Internal","(1185)", sep="\n"), "Annotated", "Yes", "790", "66.7",
  paste("Internal","(15383)", sep="\n"), "Novel", "No", "7875", "51.2",
  paste("Internal","(15383)", sep="\n"), "Novel", "Yes", "7508", "48.8",
  paste("Internal","(16568)", sep="\n"), "Total", "No", "8270", "49.9",
  paste("Internal","(16568)", sep="\n"), "Total", "Yes", "8298", "50.1",
  paste("Penultimate","(742)", sep="\n"), "Annotated", "No", "131", "17.7",
  paste("Penultimate","(742)", sep="\n"), "Annotated", "Yes", "611", "82.3",
  paste("Penultimate","(3572)", sep="\n"), "Novel", "No", "1283", "35.9",
  paste("Penultimate","(3572)", sep="\n"), "Novel", "Yes", "2289", "64.1",
  paste("Penultimate","(4314)", sep="\n"), "Total", "No", "1414", "32.8",
  paste("Penultimate","(4314)", sep="\n"), "Total", "Yes", "2900", "67.2"), byrow=T, ncol=5))
colnames(bar_data) <- c("Position", "Annotation", "Support", "Number", "Percentage")
bar_data$Percentage <- as.numeric(as.character(bar_data$Percentage))

chisq.test(matrix(as.numeric(as.character(bar_data[bar_data$Annotation=="Annotated", 4])), 2), 2)
chisq.test(matrix(as.numeric(as.character(bar_data[bar_data$Annotation=="Novel", 4])), 2), 2)
chisq.test(matrix(as.numeric(as.character(bar_data[bar_data$Annotation=="Total", 4])), 2), 2)

pdf("Human_Compare_internal_vs_penultimate.pdf")
p <- ggplot(bar_data, aes(x=Position, y=Percentage, fill=Support)) +
  geom_bar(stat = "identity", colour="black") +
  theme_bw() +
  facet_wrap(.~Annotation, scales="free") +
  scale_fill_manual(name="RNA-seq support", values=c("#999999", "#E31A1C")) +
  xlab("") +
  ylab("Percentage of exons") +
  theme(axis.text.x = element_text(size = 10, family = "Helvetica")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(text = element_text(size=18, family = "Helvetica")) +
  theme(legend.position = "top") 
print(p)
dev.off()

