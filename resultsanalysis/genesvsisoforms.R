# genes and isoforms for all patients #####
genes.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/deseq2results.csv",
                           header = T, sep = ",")
genes.all <- genes.all[,c(-1)]
colnames(genes.all)[1] <- "EnsID"
head(genes.all)

isoforms.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv",
                        header = T, sep = ",")

isoforms.all <- isoforms.all[,c(-1)]
colnames(isoforms.all)[1] <- "EnsID"
head(isoforms.all)

isoforms.all <- isoforms.all[!is.na(isoforms.all$padj),]
isoforms.sig <- isoforms.all[isoforms.all$padj < 0.05,]

genes.all$isoform.sig <- genes.all$GeneName %in% sub("-.*","",isoforms.sig$GeneName)

count <- sum(genes.all$padj < 0.05 & genes.all$isoform.sig == FALSE)

genes.all.sig <- genes.all[genes.all$padj < 0.05,]

par(mfrow = c(1,1))
boxplot(genes.all.sig$padj ~ genes.all.sig$isoform.sig, xlab = "Isoform Dysregulation", 
        ylab = "Adjusted p-value for gene dysregulation")

table(genes.all.sig$isoform.sig)



# genes and isoforms for up-responders #####

genes.up.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/up-responders/deseq2results.csv",
                        header = T, sep = ",")
genes.up.all <- genes.up.all[,c(-1)]
colnames(genes.up.all)[1] <- "EnsID"
head(genes.up.all)

isoforms.up.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv",
                           header = T, sep = ",")

isoforms.up.all <- isoforms.up.all[,c(-1)]
colnames(isoforms.up.all)[1] <- "EnsID"
head(isoforms.up.all)

isoforms.up.all <- isoforms.up.all[!is.na(isoforms.up.all$padj),]
isoforms.up.sig <- isoforms.up.all[isoforms.up.all$padj < 0.05,]

genes.up.all$isoform.sig <- genes.up.all$GeneName %in% sub("-.*","",isoforms.up.sig$GeneName)
table(genes.up.all$isoform.sig)

count <- sum(genes.up.all$padj < 0.05 & genes.up.all$isoform.sig == TRUE)

genes.up.all.sig <- genes.up.all[genes.up.all$padj < 0.05,]

par(mfrow = c(1,1))
boxplot(genes.up.all.sig$padj ~ genes.up.all.sig$isoform.sig, xlab = "Isoform Dysregulation", 
        ylab = "Adjusted p-value for gene dysregulation")

table(genes.up.all.sig$isoform.sig)



# dysregulated genes and isoforms numbers #####

genes.up.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/up-responders/deseq2results.csv",
                           header = T, sep = ",")

genes.up.all <- genes.up.all[,c(-1)]
colnames(genes.up.all)[1] <- "EnsID"
head(genes.up.all)

isoforms.up.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv",
                              header = T, sep = ",")

isoforms.up.all <- isoforms.up.all[,c(-1)]
colnames(isoforms.up.all)[1] <- "EnsID"
head(isoforms.up.all)

genes.up.all$num.isoforms <- genes.up.all$GeneName %in% sub("-.*","",isoforms.up.all$GeneName)

library(dplyr)

head(isoforms.up.all)
isoforms.up.all$GeneOnly <- sub("-.*","",isoforms.up.all$GeneName)


isoform.counts <- isoforms.up.all %>% group_by(GeneOnly) %>% summarise(counted = n())

merge <- merge(genes.up.all, isoform.counts, by.x = "GeneName", by.y = "GeneOnly", all.x = TRUE)

merge <- merge[order(merge$padj),]

merge$genesig <- merge$padj < 0.05
boxplot(merge$counted ~ merge$genesig)
table(genes.up.all$num.isoforms)

isoforms.up.all <- isoforms.up.all[!is.na(isoforms.up.all$padj),]
isoforms.up.sig <- isoforms.up.all[isoforms.up.all$padj < 0.05,]

genes.up.all$isoform.sig <- genes.up.all$GeneName %in% sub("-.*","",isoforms.up.sig$GeneName)
table(genes.up.all$isoform.sig)

count <- sum(genes.up.all$padj < 0.05 & genes.up.all$isoform.sig == TRUE)

genes.up.all.sig <- genes.up.all[genes.up.all$padj < 0.05,]

par(mfrow = c(1,1))
boxplot(genes.up.all.sig$padj ~ genes.up.all.sig$isoform.sig, xlab = "Isoform Dysregulation", 
        ylab = "Adjusted p-value for gene dysregulation")

table(genes.up.all.sig$isoform.sig)
