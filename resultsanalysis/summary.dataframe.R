# import libraries

library(tidyverse)

# import all data required to make summary data frame ######

# full list of genes

genes <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/PvR_geneCounts_all_LS_23062021.txt.txt",
                    header = T, sep = "\t")
genes <- genes[,c(1,10)]
genes$EnsIDsimp <- sub("\\..*","",genes$EnsID)
head(genes)

# check if there is just one occurrence of each gene - one of each version of each gene

gene.counts <- genes %>% group_by(EnsID) %>% summarise(counted = n())
gene.counts <- gene.counts[order(gene.counts$counted, decreasing = T),]

table(gene.counts$counted)

# jarid2 genes

jarid2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",")
jarid2 <- c(t(jarid2))

# gene pca ranking

gene.pca.rank <-  read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/resultsanalysis/topgenes.csv")
gene.pca.rank <- c(t(gene.pca.rank))

# gene DEA all patients

genes.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/deseq2results.csv",
                           header = T, sep = ",")


colnames(genes.all)[2] <- "EnsID"

# gene DEA up-responders

genes.up <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/up-responders/deseq2results.csv",
                       header = T, sep = ",")


colnames(genes.up)[2] <- "EnsID"

# gene DEA down-responders

genes.down <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/down-responders/deseq2results.csv",
                       header = T, sep = ",")

colnames(genes.down)[2] <- "EnsID"

# full list of isoforms

isoforms <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

isoforms <- isoforms[,c(1,10)]

isoforms$GeneNamesimp <- sub("-.*","",isoforms$GeneName)

# count number of isoforms per gene

isoform.counts <- isoforms %>% group_by(GeneNamesimp) %>% summarise(isoforms.full.counts = n())

isoform.counts

# isoform DEA all patients

isoforms.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv",
                         header = T, sep = ",")

isoforms.all$GeneNamesimp <- sub("-.*","",isoforms.all$GeneName)

isoforms.all.counts <- isoforms.all %>% group_by(GeneNamesimp) %>% summarise(isoforms.all.counts = n()) # how many isoforms of each gene actually tested


isoforms.all.counts.sig <- isoforms.all %>% group_by(GeneNamesimp, padj < 0.05) %>% summarise(isoforms.all.counts.sig = n())
isoforms.all.counts.sig <- isoforms.all.counts.sig[isoforms.all.counts.sig$`padj < 0.05` == TRUE,c(1,3)]

# isoform DEA up-responders

isoforms.up <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv",
                           header = T, sep = ",")

isoforms.up$GeneNamesimp <- sub("-.*","",isoforms.up$GeneName)

isoforms.up.counts <- isoforms.up %>% group_by(GeneNamesimp) %>% summarise(isoforms.up.counts = n()) # how many isoforms of each gene actually tested

isoforms.up.counts.sig <- isoforms.up %>% group_by(GeneNamesimp, padj < 0.05) %>% summarise(isoforms.up.counts.sig = n())
isoforms.up.counts.sig <- isoforms.up.counts.sig[isoforms.up.counts.sig$`padj < 0.05` == TRUE,c(1,3)]


# isoform DEA down-responders

isoforms.down <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results.csv",
                           header = T, sep = ",")

isoforms.down$GeneNamesimp <- sub("-.*","",isoforms.down$GeneName)

isoforms.down.counts <- isoforms.down %>% group_by(GeneNamesimp) %>% summarise(isoforms.down.counts = n()) # how many isoforms of each gene actually tested

isoforms.down.counts.sig <- isoforms.down %>% group_by(GeneNamesimp, padj < 0.05) %>% summarise(isoforms.down.counts.sig = n())
isoforms.down.counts.sig <- isoforms.down.counts.sig[isoforms.down.counts.sig$`padj < 0.05` == TRUE,c(1,3)]

# Compile gene section of results data frame ######

results <- genes # creating data frame for results

results$jarid2 <- results$EnsIDsimp %in% jarid2 # adding jarid2 column

results$pca.gene.rank <- match(results$EnsID, gene.pca.rank) # add gene pca rank

results$genes.all <- results$EnsID %in% genes.all$EnsID # adding whether tested in DEA of genes for all patients

results$genes.up <- results$EnsID %in% genes.up$EnsID # adding whether testing in DEA of gene for all up-responders

results$genes.down <- results$EnsID %in% genes.down$EnsID # adding whether testing in DEA of gene for all down-responders


# add in the genes.all padj

merge <- genes.all %>% 
  select(EnsID, padj, log2FoldChange)

results <- left_join(results, merge)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("padj.genes.all","LFC.genes.all")

head(results)

# add in the genes.up padj

merge <- genes.up %>% 
  select(EnsID, padj, log2FoldChange)

results <- left_join(results, merge)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("padj.genes.up","LFC.genes.up")

head(results)

# add in the genes.down padj

merge <- genes.down %>% 
  select(EnsID, padj, log2FoldChange)

results <- left_join(results, merge)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("padj.genes.down","LFC.genes.down")

head(results)

# adding in whether isoforms for gene included in full data

results$isoforms.full <- results$GeneName %in% isoforms$GeneNamesimp # adding whether there are any isoforms for the gene in question

# adding number of isoforms in full dataset, filtered all patient dataset, filtered up-responders and filtered down-responders

results <- left_join(results, isoform.counts, by = c("GeneName" = "GeneNamesimp"))
results <- left_join(results, isoforms.all.counts, by = c("GeneName" = "GeneNamesimp"))
results <- left_join(results, isoforms.up.counts, by = c("GeneName" = "GeneNamesimp"))
results <- left_join(results, isoforms.down.counts, by = c("GeneName" = "GeneNamesimp"))

# adding number of significant isoforms

results <- left_join(results, isoforms.all.counts.sig, by = c("GeneName" = "GeneNamesimp"))
results <- left_join(results, isoforms.up.counts.sig, by = c("GeneName" = "GeneNamesimp"))
results <- left_join(results, isoforms.down.counts.sig, by = c("GeneName" = "GeneNamesimp"))

head(results)

write.csv(results, "~/Documents/Semester3/Project/Results/resultsanalysis/summarydf.csv")





