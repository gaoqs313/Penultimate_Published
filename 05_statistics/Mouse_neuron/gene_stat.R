
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

##################################
##### prepare PSI data frame #####
##################################

## load tissue type
meta <- read.table("/scratch_space/qgao/Penultimate/2.data/neuron_SRR_Acc_List.txt", stringsAsFactor=F, row.names = 1)
colnames(meta) <- c("Tissue")

## load data
penul <- read.table("neuron_expr.txt", header = T, stringsAsFactor=F, row.names=1)

###########################################

Pcbp3 <- as.data.frame(t(penul[rownames(penul) == "Pcbp3", ]))
colnames(Pcbp3) <- "TPM"
Pcbp3$Stage <- meta[match(rownames(Pcbp3), rownames(meta)), 1]

Mbnl1 <- as.data.frame(t(penul[rownames(penul) == "Mbnl1", ]))
colnames(Mbnl1) <- "TPM"
Mbnl1$Stage <- meta[match(rownames(Mbnl1), rownames(meta)), 1]

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

dfp <- data_summary(Pcbp3, varname="TPM", groupnames=c("Stage"))
dfp$Stage <- factor(dfp$Stage, levels=c("DIV-8","DIV-4","DIV0","DIV1","DIV7","DIV16","DIV21","DIV28"), ordered=T)

dfm <- data_summary(Mbnl1, varname="TPM", groupnames=c("Stage"))
dfm$Stage <- factor(dfp$Stage, levels=c("DIV-8","DIV-4","DIV0","DIV1","DIV7","DIV16","DIV21","DIV28"), ordered=T)


pdf("neuron_rbp_expr.pdf", height=5)
p1 <- ggplot(dfp, aes(x=Stage, y=TPM, fill=Stage)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("#F781BF", "#999999", "#A65628", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3")) +
  geom_errorbar(aes(ymin=TPM-sd, ymax=TPM+sd), width=.2, position=position_dodge(.9)) +
  xlab("") +
  ylab("Gene expression (transcript per million)") +
  ggtitle("Pcbp3") +
  theme_bw() +
  theme(text = element_text(size=18, family = "Helvetica")) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p1)
p2 <- ggplot(dfm, aes(x=Stage, y=TPM, fill=Stage)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("#F781BF", "#999999", "#A65628", "#E41A1C", "#FF7F00", "#377EB8", "#4DAF4A", "#984EA3")) +
  geom_errorbar(aes(ymin=TPM-sd, ymax=TPM+sd), width=.2, position=position_dodge(.9)) +
  xlab("") +
  ylab("Gene expression (transcript per million)") +
  ggtitle("Mbnl1") +
  theme_bw() +
  theme(text = element_text(size=18, family = "Helvetica")) +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
print(p2)
dev.off()


