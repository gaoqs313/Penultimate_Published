
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
meta <- read.table("mouse_SRR_Acc_List.txt", stringsAsFactor=F, row.names = 1)
colnames(meta) <- c("Tissue", "Library")

## load data
internal <- read.table("internal_psi.txt", header = T, stringsAsFactor=F)

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
#13805 total
sum(internal_as$Support)
#8219 supported
table(internal_as[,c(2,4)])
#              Support
#Annotation     0    1
# Annotated   153  370
#     Novel  5433 7849

## compare penultimate with internal
bar_data <- as.data.frame(matrix(c(
  paste("Internal","(523)", sep="\n"), "Annotated", "No", "153", "29.3",
  paste("Internal","(523)", sep="\n"), "Annotated", "Yes", "370", "70.7",
  paste("Internal","(13282)", sep="\n"), "Novel", "No", "5433", "40.9",
  paste("Internal","(13282)", sep="\n"), "Novel", "Yes", "7849", "59.1",
  paste("Internal","(13805)", sep="\n"), "Total", "No", "5586", "40.5",
  paste("Internal","(13805)", sep="\n"), "Total", "Yes", "8219", "59.5",
  paste("Penultimate","(354)", sep="\n"), "Annotated", "No", "46", "13.0",
  paste("Penultimate","(354)", sep="\n"), "Annotated", "Yes", "308", "87.0",
  paste("Penultimate","(2879)", sep="\n"), "Novel", "No", "1035", "35.9",
  paste("Penultimate","(2879)", sep="\n"), "Novel", "Yes", "1844", "64.1",
  paste("Penultimate","(3233)", sep="\n"), "Total", "No", "1081", "33.4",
  paste("Penultimate","(3233)", sep="\n"), "Total", "Yes", "2152", "66.6"), byrow=T, ncol=5))
colnames(bar_data) <- c("Position", "Annotation", "Support", "Number", "Percentage")
bar_data$Percentage <- as.numeric(as.character(bar_data$Percentage))

chisq.test(matrix(as.numeric(as.character(bar_data[bar_data$Annotation=="Annotated", 4])), 2), 2)
chisq.test(matrix(as.numeric(as.character(bar_data[bar_data$Annotation=="Novel", 4])), 2), 2)
chisq.test(matrix(as.numeric(as.character(bar_data[bar_data$Annotation=="Total", 4])), 2), 2)

pdf("Compare_internal_vs_penultimate.pdf")
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

