
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
library(Rtsne)

##################################
##### prepare PSI data frame #####
##################################

## load tissue type
meta <- read.table("mouse_SRR_Acc_List.txt", stringsAsFactor=F, row.names = 1)
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

## calculate number of alternatively spliced exons
penul_as <- data.frame(Number=apply(penul_avg<90, 1, sum))
penul_as$Annotation <- penul[match(rownames(penul_as), penul$Event),2]
penul_as$Type <- penul[match(rownames(penul_as), penul$Event),3]
penul_as$Support <- ifelse(penul_as$Number > 0, 1, 0)

## output
write.table(penul_as, "mouse_support_info.txt", row.names=T, col.names=NA, quote=F, sep="\t")

## statistics
length(penul_as$Support)
#3233 total
sum(penul_as$Support)
#2152 supported
table(penul_as[,c(2,4)])
#              Support
#Annotation     0    1
# Annotated    46  308
#     Novel  1035 1844

# whether supported by RNA-seq
pie_data <- data.frame(
  state = c("without", "with"),
  amount = c(table(penul_as[,4])),
  stringsAsFactors = FALSE)
p1 <- ggplot(pie_data) +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = amount, fill = state, explode = 0), stat = 'pie') +
    scale_fill_manual(name="", values=c("#E31A1C", "#999999")) +
    coord_fixed() +
    theme_void() +
    theme(legend.position="none")

# whether annotated
pie_data1 <- data.frame(
  state = c("without", "with"),
  amount = c(table(penul_as[penul_as[,4]==1,2])),
  stringsAsFactors = FALSE)
p2 <- ggplot(pie_data1) +
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = amount, fill = state, explode = 0), stat = 'pie') +
    scale_fill_manual(name="", values=c("#c40406","#fa7374")) +
    coord_fixed() +
    theme_void() +
    theme(legend.position="none")

# stop codon position
pie_data2 <- data.frame(
  state = c(paste("Ultimate", "without RNAseq", sep='\n'),
            paste("Ultimate", "with RNAseq", sep='\n'),
            paste("Penultimate", "without RNAseq", sep='\n'),
            paste("Penultimate", "with RNAseq", sep='\n')),
  amount = c(t(table(penul_as[,3:4]))),
  stringsAsFactors = FALSE)
p3 <- ggplot(pie_data2) + 
    geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1, amount = amount, fill = state, explode = 0), stat = 'pie') + 
    scale_fill_manual(name="", values=c("#FB9A99", "grey80", "#E31A1C", "grey50")) +
    coord_fixed() +
    theme_void() +
    theme(legend.key.size = unit(1.5, 'lines'))

## generate supporting pie chart
pdf("RNAseq_support.pdf")
print(p1)
print(p2)
print(p3)
dev.off()

#################################################################
## calculate number of alternatively spliced exons in each tissue
penul_tissue <- as.data.frame(penul_avg)
penul_tissue$Type <- penul[match(rownames(penul_tissue), penul$Event),3]
penul_tissue$Type[penul_tissue$Type=="Last"] <- "Ultimate"
penul_tissue_split <- split(penul_tissue, penul_tissue$Type)
penul_tissue_split_number <- lapply(penul_tissue_split, function(x) apply(x[,1:22], 2, function(x) sum(x<90)))
penul_number_tissue <- do.call(cbind, penul_tissue_split_number)

## generate tissue supporting barplot
bar_data <- melt(penul_number_tissue)

pdf("RNAseq_tissue_support.pdf", height=6)
p <- ggplot(bar_data, aes(x=Var1, y=value, fill=Var2)) +
	geom_bar(stat = "identity", colour="black") + 
	scale_fill_manual(name = NULL, values= c("grey", "skyblue")) +
	theme_bw() +
	xlab("Tissue") +
	ylab("Number of alternatively spliced exons") +
	theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, family = "Helvetica")) +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(text = element_text(size=18, family = "Helvetica")) +
	theme(legend.position=c(0.15, 0.9))
print(p)
dev.off()

###############################################################
## calculate number of AS exons vs number of supporting tissues
penul_tissue <- as.data.frame(penul_avg)

## get expressed in all and alternative in at least 1 tissue
penul_tissue_exp <- penul_tissue[apply(penul_tissue<1000, 1, all) & apply(penul_tissue<90, 1, any), ]

## calculate number of supporting tissues
penul_tissue_exp_num <- apply(penul_tissue_exp<90, 1, sum)

## generate tissue supporting histogram
hist_data <- as.data.frame(table(penul_tissue_exp_num))

pdf("RNAseq_tissue_distribution.pdf", height=6)
p <- ggplot(hist_data, aes(x=penul_tissue_exp_num, y=Freq)) +
        geom_bar(stat = "identity", colour="black", fill="grey") +
        theme_bw() +
        xlab("Number of supporting tissues") +
        ylab("Number of alternatively spliced exons") +
        theme(axis.text.x = element_text(size = 9, family = "Helvetica")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        theme(text = element_text(size=18, family = "Helvetica"))
print(p)
dev.off()

##############################################
##### generate tissue-specificity figure #####
##############################################

## generate tissue specificity
psi <- column_to_rownames(psi)
psis <- psi[rownames(psi) %in% rownames(penul_as[penul_as[,4]==1,]), ]
#2152

## expressed in at least 50% (>=55) of samples
psie <- psis[apply(!is.na(psis), 1, sum) >=55, ]
#1224

## batch effect correction
modcombat = model.matrix(~1, data=meta)
psib = ComBat(dat=as.matrix(psie), batch=meta$Library, mod=modcombat, par.prior=TRUE) %>% as.data.frame()

## calculate mad
psis_mad <- apply(psib, 1, function(x) mad(x, na.rm=T))

## calculate pairwise weighted correlation
psis_wcor <- wtd.cors(psib, weight=psis_mad)
myDist <- as.dist(1 - psis_wcor)

## tSNE
set.seed(0)
psis_tsne = Rtsne(myDist, dims = 2, perplexity = 5, theta = 0, is_distance = T, pca = F, max_iter = 5000, check_duplicates = F ) # Run TSNE

## prepare color
tissueCol=c()
{
  tissueCol["Adrenal"]="gold2"
  tissueCol["Bladder"]="pink"
  tissueCol["Cerebellum"]="#1F78B5"
  tissueCol["Colon"]="purple"
  tissueCol["Cortex"]='red4'
  tissueCol["Duodenum"]="darkgoldenrod4"
  tissueCol["FrontalLobe"]="#A8DD00"
  tissueCol["GenitalFatPad"]="black"
  tissueCol["Heart"]="#CCCC33"
  tissueCol["Kidney"]='#A65628'
  tissueCol["LgIntestine"]= "skyblue"
  tissueCol["Liver"]="#1E90FF"
  tissueCol["Lung"]="#3E9F32"
  tissueCol["MammaryGland"]="magenta3"
  tissueCol["Ovary"]="#1E90FF"
  tissueCol["Placenta"]='blue3'
  tissueCol["SmIntestine"]="grey40"
  tissueCol["Spleen"]="lightslateblue"
  tissueCol["Stomach"]="seagreen2"
  tissueCol["SubcFatPad"]="orangered"
  tissueCol["Testis"]="#66C2A6"
  tissueCol["Thymus"]='cyan'
}

## plot tSNE
tsneDF = meta %>% mutate(X=psis_tsne$Y[,1], Y=psis_tsne$Y[,2]) 

pdf("tSNE_RNAseq.pdf")
p <- ggplot() + xlab("tSNE-1") + ylab("tSNE-2")+
    geom_point(data=tsneDF, aes(X, Y, color=Tissue), size=1)+
    theme_bw() +
    scale_color_manual(values = tissueCol) +
    theme(axis.text = element_blank(), axis.ticks=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(text = element_text(size=20, family = "Helvetica"))
centers <- tsneDF %>% group_by(Tissue) %>% summarise(x=mean(x=X), y=mean(x=Y))
p <- p + geom_point(data=centers, mapping=aes(x=x, y=y), size=0, alpha=0) + geom_text(data=centers, mapping=aes(x=x, y=y,label=Tissue), size=5, family="Helvetica", color="black") + theme(legend.position = "none")
print(p)
dev.off()

