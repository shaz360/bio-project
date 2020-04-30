#set directory for RBBP5 vs N2 to dropbox location
setwd("/Users/Shazaib/Dropbox/UG BioInfo 2019-20/SHAZT1/R5 vs N2")
#check if new directory is set
getwd()

#READ the R5VSN2 csv file and label R5VSN2 
dat_R5VSN2<-read.csv("./R5VSN2.csv", header = TRUE, stringsAsFactors = FALSE)

#CHECK the csv to see if all the data is there
# - head and tail = see first 6 and last 6 rows
head(dat_R5VSN2)
tail(dat_R5VSN2)
# - see names and number of rows and their class
names(dat_R5VSN2)
nrow (dat_R5VSN2)
sapply(dat_R5VSN2, class)

################FC2/PVAL0.01##########################

#SUBSET the differentially expressed (DE) genes
# by Fold change 2 and p value 0.01
dat_R5VSN2UPDOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=1 & dat_R5VSN2$pval<=0.01|dat_R5VSN2$log2FoldChange<=-1 & dat_R5VSN2$pval<=0.01, )
#CHECK
head(dat_R5VSN2UPDOWN)
tail(dat_R5VSN2UPDOWN)
names(dat_R5VSN2UPDOWN)
nrow(dat_R5VSN2UPDOWN)
sapply(dat_R5VSN2UPDOWN, class)

#SUBSET only upregulated genes
dat_R5VSN2UP<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=1 & dat_R5VSN2$pval<=0.01)
#CHECK
head(dat_R5VSN2UP)
tail(dat_R5VSN2UP)
names(dat_R5VSN2UP)
nrow(dat_R5VSN2UP)
sapply(dat_R5VSN2UP, class)

#SUBSET only downregulated genes
dat_R5VSN2DOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange<=-1 & dat_R5VSN2$pval<=0.01)
#CHECK
head(dat_R5VSN2DOWN)
tail(dat_R5VSN2DOWN)
names(dat_R5VSN2DOWN)
nrow(dat_R5VSN2DOWN)
sapply(dat_R5VSN2DOWN, class)

#Create 3 separate csv files containing the subset data
write.table(dat_R5VSN2UPDOWN, file="YAR5VSN2UPDOWN.csv", row.names = F, sep=",")
write.table(dat_R5VSN2UP, file="YAR5VSN2UP.csv", row.names = F, sep=",")
write.table(dat_R5VSN2DOWN, file="YAR5VSN2DOWN.csv", row.names = F, sep=",")

##################FC 1.5/PVAL0.05#################################

#SUBSET the differentially expressed (DE) genes
# by Fold change 1.5 and p value 0.05 
dat_1p50p05YAR5VSN2UPDOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=0.58496250072 & dat_R5VSN2$pval<=0.05|dat_R5VSN2$log2FoldChange<=-0.58496250072 & dat_R5VSN2$pval<=0.05, )

#CHECK
head(dat_1p50p05YAR5VSN2UPDOWN)
tail(dat_1p50p05YAR5VSN2UPDOWN)
names(dat_1p50p05YAR5VSN2UPDOWN)
nrow(dat_1p50p05YAR5VSN2UPDOWN)
sapply(dat_1p50p05YAR5VSN2UPDOWN, class)

#SUBSET only upregulated genes
dat_1p50p05YAR5VSN2UP<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=0.58496250072 & dat_R5VSN2$pval<=0.05)
#check
head(dat_1p50p05YAR5VSN2UP)
tail(dat_1p50p05YAR5VSN2UP)
names(dat_1p50p05YAR5VSN2UP)
nrow(dat_1p50p05YAR5VSN2UP)
sapply(dat_1p50p05YAR5VSN2UP, class)

#SUBSET only downregulated genes
dat_1p50p05YAR5VSN2DOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange<=-0.58496250072 & dat_R5VSN2$pval<=0.05)
#CHECK
head(dat_1p50p05YAR5VSN2DOWN)
tail(dat_1p50p05YAR5VSN2DOWN)
names(dat_1p50p05YAR5VSN2DOWN)
nrow(dat_1p50p05YAR5VSN2DOWN)
sapply(dat_1p50p05YAR5VSN2DOWN, class)

###############FC1.5/PVAL0.01###########

#SUBSET the differentially expressed (DE) genes
# by Fold change 1.5 and p value 0.01 
dat_1p50p01YAR5VSN2UPDOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=0.58496250072 & dat_R5VSN2$pval<=0.01|dat_R5VSN2$log2FoldChange<=-0.58496250072 & dat_R5VSN2$pval<=0.01, )
#CHECK
head(dat_1p50p01YAR5VSN2UPDOWN)
tail(dat_1p50p01YAR5VSN2UPDOWN)
names(dat_1p50p01YAR5VSN2UPDOWN)
nrow(dat_1p50p01YAR5VSN2UPDOWN)
sapply(dat_1p50p01YAR5VSN2UPDOWN, class)

#SUBSET only upregulated genes
dat_1p50p01YAR5VSN2UP<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=0.58496250072 & dat_R5VSN2$pval<=0.01)
#check
head(dat_1p50p01YAR5VSN2UP)
tail(dat_1p50p01YAR5VSN2UP)
names(dat_1p50p01YAR5VSN2UP)
nrow(dat_1p50p01YAR5VSN2UP)
sapply(dat_1p50p01YAR5VSN2UP, class)

#SUBSET only downregulated genes
dat_1p50p01YAR5VSN2DOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange<=-0.58496250072 & dat_R5VSN2$pval<=0.01)
#CHECK
head(dat_1p50p01YAR5VSN2DOWN)
tail(dat_1p50p01YAR5VSN2DOWN)
names(dat_1p50p01YAR5VSN2DOWN)
nrow(dat_1p50p01YAR5VSN2DOWN)
sapply(dat_1p50p01YAR5VSN2DOWN, class)

##############FC2/PVAL0.05###############

#SUBSET the differentially expressed (DE) genes
# by Fold change 2 and p value<0.05 
dat_20p05YAR5VSN2UPDOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=1 & dat_R5VSN2$pval<=0.05|dat_R5VSN2$log2FoldChange<=-1 & dat_R5VSN2$pval<=0.05, )
#CHECK
head(dat_20p05YAR5VSN2UPDOWN)
tail(dat_20p05YAR5VSN2UPDOWN)
names(dat_20p05YAR5VSN2UPDOWN)
nrow(dat_20p05YAR5VSN2UPDOWN)
sapply(dat_20p05YAR5VSN2UPDOWN, class)

#SUBSET only upregulated genes
dat_20p05YAR5VSN2UP<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=1 & dat_R5VSN2$pval<=0.05)
#check
head(dat_20p05YAR5VSN2UP)
tail(dat_20p05YAR5VSN2UP)
names(dat_20p05YAR5VSN2UP)
nrow(dat_20p05YAR5VSN2UP)
sapply(dat_20p05YAR5VSN2UP, class)

#SUBSET only downregulated genes
dat_20p05YAR5VSN2DOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange<=-1 & dat_R5VSN2$pval<=0.05)
#CHECK
head(dat_20p05YAR5VSN2DOWN)
tail(dat_20p05YAR5VSN2DOWN)
names(dat_20p05YAR5VSN2DOWN)
nrow(dat_20p05YAR5VSN2DOWN)
sapply(dat_20p05YAR5VSN2DOWN, class)

#########FC3/PVAL0.001################

#SUBSET the differentially expressed (DE) genes
# by Fold change 3 and p value 0.001 
#label the results DiffExpYAR5VSN2
dat_30p001YAR5VSN2UPDOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=1.58496250072 & dat_R5VSN2$pval<=0.001|dat_R5VSN2$log2FoldChange<=-1.58496250072 & dat_R5VSN2$pval<=0.001, )
#CHECK
head(dat_30P001YAR5VSN2UPDOWN)
tail(dat_30P001YAR5VSN2UPDOWN)
names(dat_30P001YAR5VSN2UPDOWN)
nrow(dat_30P001YAR5VSN2UPDOWN)
sapply(dat_30P001YAR5VSN2UPDOWN, class)

#SUBSET only upregulated genes
dat_30P001YAR5VSN2UP<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange>=1.58496250072 & dat_R5VSN2$pval<=0.001)
#check
head(dat_30P001YAR5VSN2UP)
tail(dat_30P001YAR5VSN2UP)
names(dat_30P001YAR5VSN2UP)
nrow(dat_30P001YAR5VSN2UP)
sapply(dat_30P001YAR5VSN2UP, class)

#SUBSET only downregulated genes
dat_30P001YAR5VSN2DOWN<-subset(dat_R5VSN2, dat_R5VSN2$log2FoldChange<=-1.58496250072 & dat_R5VSN2$pval<=0.001)
#CHECK
head(dat_30P001YAR5VSN2DOWN)
tail(dat_30P001YAR5VSN2DOWN)
names(dat_30P001YAR5VSN2DOWN)
nrow(dat_30P001YAR5VSN2DOWN)
sapply(dat_30P001YAR5VSN2DOWN, class)

###################################################################

#VOLCANO PLOT
#Turn on enhancedvolcano package
library(EnhancedVolcano)

#create volc plot
EnhancedVolcano(dat_R5VSN2,
                lab = rownames(dat_R5VSN2), #required
                x = "log2FoldChange", #r
                y = "pval", #r
                xlim = c(-4, 4),
                ylim = c(0,6),
                FCcutoff = 1,
                pCutoff = 10e-2,
                pointSize = 1,
                labSize = 2.5, 
                selectLab = 0,
                legendVisible = FALSE,
                title = "EMBRYO rbbp-5(-) vs wild type (N2)",
                subtitle = "p-value against fold change",
                caption= paste0('Total= ', nrow(dat_R5VSN2), '  Misregulated= ', nrow(dat_R5VSN2UPDOWN), '  Upregulated= ', nrow(dat_R5VSN2UP), '  Downregulated= ', nrow(dat_R5VSN2DOWN)), 
                titleLabSize = 14,
                axisLabSize = 12,
                legend = 1,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 1)

