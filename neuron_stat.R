
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
meta <- read.table("neuron_SRR_Acc_List.txt", stringsAsFactor=F, row.names = 1)
colnames(meta) <- c("Stage")

## load data
psi <- read.table("neuron_psi.txt", header = T, stringsAsFactor=F, row.names = 1)

#########################################################
##### generate development stage-specificity figure #####
#########################################################

## expressed in at least 1 of samples
psie <- psi[apply(!is.na(psi), 1, sum) > 0, ]
#2458

## calculate mad
psi_mad <- apply(psie, 1, function(x) mad(x, na.rm=T))

## calculate pairwise weighted correlation
psi_wcor <- wtd.cors(psie, weight=psi_mad)
myDist <- as.dist(1 - psi_wcor)

## tSNE
set.seed(0)
psi_tsne = Rtsne(myDist, dims = 2, perplexity = 4, theta = 0, is_distance = T, pca = F, max_iter = 5000, check_duplicates = F ) # Run TSNE

## plot tSNE
tsneDF = meta %>% mutate(X=psi_tsne$Y[,1], Y=psi_tsne$Y[,2]) 
tsneDF$Stage <- factor(tsneDF$Stage, levels=c("DIV-8","DIV-4","DIV0","DIV1","DIV7","DIV16","DIV21","DIV28"), ordered=T)

pdf("tSNE_neuron.pdf")
p <- ggplot() + xlab("tSNE-1") + ylab("tSNE-2")+
    geom_point(data=tsneDF, aes(X, Y, color=Stage), size=2)+
    theme_bw() +
    scale_color_manual(values = c("#F781BF", "#999999", "#A65628", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3")) +
    xlim(min(tsneDF$X)-10, max(tsneDF$X)+10) +
    ylim(min(tsneDF$Y)-10, max(tsneDF$Y)+10) +
    theme(axis.text = element_blank(), axis.ticks=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(text = element_text(size=20, family = "Helvetica"))
centers <- tsneDF %>% group_by(Stage) %>% summarise(x=mean(x=X), y=mean(x=Y))
p <- p + geom_point(data=centers, mapping=aes(x=x, y=y), size=0, alpha=0) + geom_text(data=centers, mapping=aes(x=x, y=y,label=Stage), size=5, family="Helvetica", color="black") + theme(legend.position = "none")
print(p)
dev.off()

##################################################
##### find extreme examples for illustration #####
##################################################

## melt
psie$Event <- rownames(psie)
psi_melt <- melt(psie)

## add stage
psi_melt$Stage <- meta[match(psi_melt$variable, rownames(meta)), 1]

## split by stage
psi_split <- split(psi_melt, psi_melt$Stage)

## calculate average psi for each tissue
psi_avg_list = list()
for(i in 1:length(psi_split))
{
  tmp = split(psi_split[[i]], psi_split[[i]][,1])
  psi_avg_list[[i]] = lapply(tmp, function(x) mean(x$value, na.rm=T))
}
psi_avg <- do.call(cbind, psi_avg_list)
colnames(psi_avg) <- sort(unique(meta$Stage))
psi_avg[psi_avg=="NaN"] <- 1000

## express in all stages
psi_all <- psi_avg[apply(psi_avg<1000, 1, all), ]

## reorder
psi_order <- psi_all[, c(2,1,3,4,8,5:7)]
psi_order <- apply(psi_order, 2, as.numeric)
rownames(psi_order) <- rownames(psi_all)

## increase
psi_increase <- psi_order[apply(psi_order[,2:8]-psi_order[,1:7]>0, 1, all),, drop=F]
psi_increase <- psi_order[apply(psi_order[,3:8]-psi_order[,2:7]>0, 1, all),, drop=F]

## decrease
psi_decrease <- psi_order[apply(psi_order[,2:8]-psi_order[,1:7]<0, 1, all),, drop=F]

## get example
examp_inc <- as.data.frame(t(psi[rownames(psi) == "14:121545966:121546104:-@14:121544576:121544583:-@14:121542046:121543625:-",]))
colnames(examp_inc) <- "PSI"
examp_inc$Stage <- meta[match(rownames(examp_inc), rownames(meta)), 1] 
examp_dec <- as.data.frame(t(psi[rownames(psi) == "11:95046921:95047046:-@11:95046449:95046590:-@11:95044474:95045763:-",]))
colnames(examp_dec) <- "PSI"
examp_dec$Stage <- meta[match(rownames(examp_dec), rownames(meta)), 1]

## data summary
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
 return(data_sum)
}

df_inc <- data_summary(examp_inc, varname="PSI", groupnames=c("Stage"))
df_dec <- data_summary(examp_dec, varname="PSI", groupnames=c("Stage"))
df_inc$Stage <- factor(df_inc$Stage, levels=c("DIV-8","DIV-4","DIV0","DIV1","DIV7","DIV16","DIV21","DIV28"), ordered=T)
df_dec$Stage <- factor(df_dec$Stage, levels=c("DIV-8","DIV-4","DIV0","DIV1","DIV7","DIV16","DIV21","DIV28"), ordered=T)

pdf("neuron_example.pdf", height=5)
p1 <- ggplot(df_inc, aes(x=Stage, y=PSI, fill=Stage)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("#F781BF", "#999999", "#A65628", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3")) +
  geom_errorbar(aes(ymin=PSI-sd, ymax=PSI+sd), width=.2, position=position_dodge(.9)) +
  xlab("") +
  ylab("Percent Splice In") +
  ggtitle("Dock9 (ENSMUST00000212376)") +
  theme_bw() +
  theme(text = element_text(size=18, family = "Helvetica")) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p1)
p2 <- ggplot(df_dec, aes(x=Stage, y=PSI, fill=Stage)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("#F781BF", "#999999", "#A65628", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3")) +
  geom_errorbar(aes(ymin=PSI-sd, ymax=PSI+sd), width=.2, position=position_dodge(.9)) +
  xlab("") +
  ylab("Percent Splice In") +
  ggtitle("Itga3 (ENSMUST00000107739)") +
  theme_bw() +
  theme(text = element_text(size=18, family = "Helvetica")) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()

