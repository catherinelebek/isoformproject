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

gene.counts <- genes %>% group_by(EnsIDsimp) %>% summarise(counted = n())
gene.counts <- gene.counts[order(gene.counts$counted, decreasing = T),]

length(unique(genes$EnsIDsimp))

# only duplicates belong to genes within a pseudo-autosomal region (PAR_Y)

dupes <- gene.counts[gene.counts$counted == 2,1]
dupes
nrow(genes[grep("PAR_Y",genes$EnsID),])

genes[genes$EnsIDsimp %in% dupes$EnsIDsimp,]

# jarid2 genes

jarid2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",")
jarid2 <- c(t(jarid2))

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

# add in mapping data

mapping <- read.delim("~/Documents/Semester3/Project/InputData/output.txt", header = F, sep = "\t") # gene-to-transcript IDs

colnames(mapping) <- c("Gene.EnsID","Transcript.EnsID")

mapping$Gene.EnsID.Simp <- sub("\\..*","",mapping$Gene.EnsID)
mapping$Transcript.EnsID.Simp <- sub("\\..*","",mapping$Transcript.EnsID)

# full list of isoforms

isoforms <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

isoforms <- isoforms[,c(1,10)]

isoforms$GeneNamesimp <- sub("-.*","",isoforms$GeneName)
isoforms$EnsIDsimp <- sub("\\..*","",isoforms$EnsID)

head(isoforms)
head(mapping)

isoforms <- left_join(isoforms, mapping, by = c("EnsID" = "Transcript.EnsID"))

length(unique(isoforms$EnsID))
isoform.counts <- isoforms %>% group_by(EnsIDsimp) %>% summarise(counted = n())
isoform.counts <- isoform.counts[order(isoform.counts$counted, decreasing = T),]
head(isoform.counts)

# only duplicates belong to isoforms from genes within a pseudo-autosomal region (PAR_Y)

dupes <- isoform.counts[isoform.counts$counted == 2,1]
dupes
nrow(isoforms[grep("PAR_Y",isoforms$EnsID),])

isoforms[isoforms$EnsIDsimp %in% dupes$EnsIDsimp,]

# isoform DEA all patients

isoforms.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv",
                         header = T, sep = ",")

isoforms.all$GeneNamesimp <- sub("-.*","",isoforms.all$GeneName)

head(isoforms.all)

isoforms.all.counts <- isoforms.all %>% group_by(GeneNamesimp) %>% summarise(isoforms.all.counts = n()) # how many isoforms of each gene actually tested
isoforms.all.counts.sig <- isoforms.all %>% group_by(GeneNamesimp, padj < 0.05, abs(log2FoldChange) > 1) %>% summarise(isoforms.all.counts.sig = n())
isoforms.all.counts.sig <- isoforms.all.counts.sig[isoforms.all.counts.sig$`padj < 0.05` == TRUE & 
                                                   isoforms.all.counts.sig$`abs(log2FoldChange) > 1` == TRUE,
                                                    c(1,4)]

# isoform DEA up-responders

isoforms.up <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv",
                           header = T, sep = ",")

isoforms.up$GeneNamesimp <- sub("-.*","",isoforms.up$GeneName)

isoforms.up.counts <- isoforms.up %>% group_by(GeneNamesimp) %>% summarise(isoforms.up.counts = n()) # how many isoforms of each gene actually tested

isoforms.up.counts.sig <- isoforms.up %>% group_by(GeneNamesimp, padj < 0.05, abs(log2FoldChange) > 1) %>% summarise(isoforms.up.counts.sig = n())
isoforms.up.counts.sig <- isoforms.up.counts.sig[isoforms.up.counts.sig$`padj < 0.05` == TRUE & 
                                                 isoforms.up.counts.sig$`abs(log2FoldChange) > 1` == TRUE,
                                                 c(1,4)]

# isoform DEA down-responders

isoforms.down <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results.csv",
                           header = T, sep = ",")

isoforms.down$GeneNamesimp <- sub("-.*","",isoforms.down$GeneName)

isoforms.down.counts <- isoforms.down %>% group_by(GeneNamesimp) %>% summarise(isoforms.down.counts = n()) # how many isoforms of each gene actually tested

isoforms.down.counts.sig <- isoforms.down %>% group_by(GeneNamesimp, padj < 0.05, abs(log2FoldChange) > 1) %>% summarise(isoforms.down.counts.sig = n())
isoforms.down.counts.sig <- isoforms.down.counts.sig[isoforms.down.counts.sig$`padj < 0.05` == TRUE &
                                                     isoforms.down.counts.sig$`abs(log2FoldChange) > 1` == TRUE,
                                                     c(1,4)]

# Compile gene section of results data frame ######

results <- genes # creating data frame for results

results$jarid2 <- results$EnsIDsimp %in% jarid2 # adding jarid2 column

results$genes.all <- results$EnsID %in% genes.all$EnsID # adding whether tested in DEA of genes for all patients

results$genes.up <- results$EnsID %in% genes.up$EnsID # adding whether testing in DEA of gene for all up-responders

results$genes.down <- results$EnsID %in% genes.down$EnsID # adding whether testing in DEA of gene for all down-responders


# add in the genes.all padj

merge <- genes.all %>% dplyr::select(EnsID, padj, log2FoldChange)

results <- left_join(results, merge)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("padj.genes.all","LFC.genes.all")

head(results)

# add in the genes.up padj

merge <- genes.up %>% dplyr::select(EnsID, padj, log2FoldChange)

results <- left_join(results, merge)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("padj.genes.up","LFC.genes.up")

head(results)

# add in the genes.down padj

merge <- genes.down %>% dplyr::select(EnsID, padj, log2FoldChange)

results <- left_join(results, merge)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("padj.genes.down","LFC.genes.down")

head(results)

 # adding in whether isoforms for gene included in full data

test <- isoforms[!is.na(isoforms$Gene.EnsID),]
table(test$Gene.EnsID %in% results$EnsID)
table(results$EnsID %in% test$Gene.EnsID)

results$isoforms.full <- results$EnsID %in% isoforms$Gene.EnsID # adding whether there are any isoforms for the gene in question

test1 <- isoforms[!isoforms$Gene.EnsID %in% mapping$Gene.EnsID,3]
test2 <- results[results$isoforms.full == FALSE,2]

head(isoforms)
head(results)

results$isoforms.full <- ifelse(results$isoforms.full == FALSE,
                                      results$GeneName %in% isoforms$GeneNamesimp,
                                      results$isoforms.full)

table(results$isoforms.full)

# adding number of isoforms in full dataset, filtered all patient dataset, filtered up-responders and filtered down-responders

head(isoform.counts)
head(results)


results <- full_join(results, isoform.counts, by = c("GeneName" = "GeneNamesimp"))
results <- full_join(results, isoforms.all.counts, by = c("GeneName" = "GeneNamesimp"))
results <- full_join(results, isoforms.up.counts, by = c("GeneName" = "GeneNamesimp"))
results <- full_join(results, isoforms.down.counts, by = c("GeneName" = "GeneNamesimp"))


# adding number of significant isoforms

results <- left_join(results, isoforms.all.counts.sig, by = c("GeneName" = "GeneNamesimp"))
results <- left_join(results, isoforms.up.counts.sig, by = c("GeneName" = "GeneNamesimp"))
results <- left_join(results, isoforms.down.counts.sig, by = c("GeneName" = "GeneNamesimp"))

head(results)

write.csv(results, "~/Documents/Semester3/Project/Results/resultsanalysis/summarydf.csv")

# pulling out the stats I need to compare isoform and gene DE #####

# genes significant

idx <- !is.na(results$padj.genes.all)
res <- results[idx,]
unique(length(res$GeneName))

# all patients

idx <- !is.na(results$padj.genes.all)
res <- results[idx,]

unique(length(res$GeneName))

table(res$padj.genes.all < 0.05 & res$jarid2 == FALSE)
table(res$padj.genes.all < 0.05 & res$LFC.genes.all < 0 & res$jarid2 == FALSE)
table(res$padj.genes.all < 0.05 & res$LFC.genes.all < 0 & res$jarid2 == TRUE)
table(res$padj.genes.all < 0.05 & abs(res$LFC.genes.all) > 1 & res$jarid2 == FALSE)

unique(length(res$GeneName))

test <- res[is.na(res$isoforms.all.counts.sig) & !is.na(res$isoforms.all.counts),]
length(unique(test$GeneName))
test <- res[is.na(res$isoforms.all.counts.sig) & is.na(res$isoforms.all.counts),]
length(unique(test$GeneName))

test <- res[!is.na(res$isoforms.all.counts.sig),]
length(unique(test$GeneName))

# isoforms significant

idx <- !is.na(results$isoforms.all.counts.sig)
res <- results[idx,]

length(unique(res$GeneName))

# genes tested

idx <- !is.na(res$padj.genes.all) & res$genes.all == TRUE
res1 <- res[idx,]

length(unique(res1$GeneName))

# genes tested and significant

res2 <- res1[res1$padj.genes.all < 0.05 & abs(res1$LFC.genes.all) > 1,]

length(unique(res2$GeneName))

# genes tested and not significant

res3 <- res1[res1$padj.genes.all > 0.05 | abs(res1$LFC.genes.all) < 1,]

length(unique(res3$GeneName))

sum(res3$GeneName %in% res2$GeneName)

# genes not tested

idx <- is.na(res$padj.genes.all) | res$padj.genes.all == FALSE
res4 <- res[idx,]
res4$check <- res4$GeneName %in% res1$GeneName
table(res4$check)

sum(res4$GeneName %in% res2$GeneName)

head(res2)

genes.all[genes.all$GeneName == "CPEB1",]
results[results$GeneName == "CPEB1",]


# not significant in either but still tested

idx <- !is.na(results$genes.all) & !is.na(results$isoforms.all.counts) & (results$padj.genes.all > 0.05 |
       abs(results$LFC.genes.all) < 1) & is.na(results$isoforms.all.counts.sig)

res <- results[idx,]

length(unique(res$GeneName))
