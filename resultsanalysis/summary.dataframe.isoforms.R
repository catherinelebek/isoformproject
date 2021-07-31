# import libraries

library(tidyverse)
library(VennDiagram)
library(ggvenn)



# full list of isoforms

isoforms <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

isoforms <- isoforms[,c(1,10)]

# add in mapping data

mapping <- read.delim("~/Documents/Semester3/Project/InputData/output.txt", header = F, sep = "\t") # gene-to-transcript IDs

colnames(mapping) <- c("Gene.EnsID","Transcript.EnsID")

mapping$Gene.EnsID.Simp <- sub("\\..*","",mapping$Gene.EnsID)
mapping$Transcript.EnsID.Simp <- sub("\\..*","",mapping$Transcript.EnsID)

table(isoforms$EnsID %in% mapping$Transcript.EnsID)

results <- left_join(isoforms, mapping, by = c("EnsID" = "Transcript.EnsID"))

# isoform results all patients

isoforms.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results_merge.csv",
                           header = T, sep = ",")

# isoform results up-responders

isoforms.up <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results_merge.csv",
                           header = T, sep = ",")

# isoform results down-responders

isoforms.down <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results_merge.csv",
                           header = T, sep = ",")

# check that the isoforms in the DEA results appear in the full list

all(isoforms.all$Row.names %in% isoforms$EnsID)
all(isoforms.up$Row.names %in% isoforms$EnsID)
all(isoforms.down$Row.names %in% isoforms$EnsID)

# results data frame

results$tested.all <- isoforms$EnsID %in% isoforms.all$Row.names
results$tested.up <- isoforms$EnsID %in% isoforms.up$Row.names
results$tested.down <- isoforms$EnsID %in% isoforms.down$Row.names

head(results)

colnames(results)[1] <- "Transcript.EnsID"

# read in JARID2 data

jarid2.gene <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                          header = F, sep = ",") # JARID2 gene IDs

jarid2.gene <- c(t(jarid2.gene))

jarid2.tss <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/jarid.tss.transcripts.csv",
                       header = T) # JARID2 TSS transcript IDs

jarid2.tss <- c(t(jarid2.tss))
jarid2.tss <- sub("\\..*","",jarid2.tss)

# add JARID2 data to results data frame

results$jarid2.gene <- results$Gene.EnsID.Simp %in% jarid2.gene

results$jarid2.tss <- results$Transcript.EnsID.Simp %in% jarid2.tss

head(results)

# make sure all the isoforms are unique

nrow(results)
length(unique(results$Transcript.EnsID))

nrow(isoforms.all)
length(unique(isoforms.all$Row.names))

nrow(isoforms.up)
length(unique(isoforms.up$Row.names))

nrow(isoforms.down)
length(unique(isoforms.down$Row.names))

# merge results

# all patients

results <- merge(results, isoforms.all[,c("Row.names","log2FoldChange","padj")], by.x = "Transcript.EnsID", by.y = "Row.names", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.all","padj.all")

# up-responders

results <- merge(results, isoforms.up[,c("Row.names","log2FoldChange","padj")], by.x = "Transcript.EnsID", by.y = "Row.names", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.up","padj.up")

# down-responders

results <- merge(results, isoforms.down[,c("Row.names","log2FoldChange","padj")], by.x = "Transcript.EnsID", by.y = "Row.names", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.down","padj.down")

head(results)

# import the gene-level data

# all patients

genes.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/deseq2results.csv",
                        header = T, sep = ",")

colnames(genes.all)[2] <- "EnsID"

# up-responders

genes.up <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/up-responders/deseq2results.csv",
                       header = T, sep = ",")

colnames(genes.up)[2] <- "EnsID"

# down-responders

genes.down <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/down-responders/deseq2results.csv",
                         header = T, sep = ",")

colnames(genes.down)[2] <- "EnsID"

# make sure all the genes in the gene dea data are unique

nrow(genes.all)
length(unique(genes.all$EnsID))

nrow(genes.up)
length(unique(genes.up$EnsID))

nrow(genes.down)
length(unique(genes.down$EnsID))

# merge results

results <- merge(results, genes.all[,c("EnsID","log2FoldChange","padj")], by.x = "Gene.EnsID", by.y = "EnsID", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.all.genes","padj.all.genes")


results <- merge(results, genes.up[,c("EnsID","log2FoldChange","padj")], by.x = "Gene.EnsID", by.y = "EnsID", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.up.genes","padj.up.genes")


results <- merge(results, genes.down[,c("EnsID","log2FoldChange","padj")], by.x = "Gene.EnsID", by.y = "EnsID", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.down.genes","padj.down.genes")

# rearrange columns

results <- results[,c(2,8,1,7,3,9,10,4,5,6,11:22)]
head(results)

# plot a venn diagram

# isoforms tested

ggplot(results, aes(A = tested.all, B = tested.up, C = tested.down)) +
  geom_venn()

# jarid2 genes vs. jarid2 genes & tss

ggplot(results, aes(A = jarid2.gene, B = jarid2.tss)) +
  geom_venn()

# significant results

results$sig.all <- ifelse(is.na(results$padj.all), FALSE, ifelse(results$padj.all < 0.05, TRUE, FALSE))
results$sig.up <- ifelse(is.na(results$padj.up), FALSE, ifelse(results$padj.up < 0.05, TRUE, FALSE))
results$sig.down <- ifelse(is.na(results$padj.down), FALSE, ifelse(results$padj.down < 0.05, TRUE, FALSE))

ggplot(results, aes(A = sig.all, B = sig.up, C = sig.down)) +
  geom_venn()

# jarid2 significant results

ggplot(results, aes(A = jarid2.tss, B = sig.all, C = jarid2.gene)) +
  geom_venn()

ggplot(results, aes(A = jarid2.tss, B = sig.up, C = jarid2.gene)) +
  geom_venn()

ggplot(results, aes(A = jarid2.tss, B = sig.down, C = jarid2.gene)) +
  geom_venn()

# genes also significant

results$sig.all.gene <- ifelse(is.na(results$padj.all.genes), FALSE, ifelse(results$padj.all.genes < 0.05, TRUE, FALSE))
results$sig.up.gene <- ifelse(is.na(results$padj.up.genes), FALSE, ifelse(results$padj.up.genes < 0.05, TRUE, FALSE))
results$sig.down.gene <- ifelse(is.na(results$padj.down.genes), FALSE, ifelse(results$padj.down.genes < 0.05, TRUE, FALSE))

ggplot(results, aes(A = sig.all, B = sig.all.gene)) +
  geom_venn()

ggplot(results, aes(A = sig.up, B = sig.up.gene)) +
  geom_venn()

ggplot(results, aes(A = sig.down, B = sig.down.gene)) +
  geom_venn()


# summarise at the gene level

results$threshold.all <- results$sig.all == TRUE & abs(results$LFC.all) > 1
results$threshold.up <- results$sig.up == TRUE & abs(results$LFC.up) > 1
results$threshold.down <- results$sig.down == TRUE & abs(results$LFC.down) > 1

results$match <- results$threshold.up != results$threshold.down

sum(is.na(results$threshold.down))

gene.counts <- results %>% group_by(Gene.EnsID) %>% summarise(total.isoforms = n(), threshold.all = sum(threshold.all), 
                                                              threshold.up = sum(threshold.up), threshold.down = sum(threshold.down),
                                                              sig.all.gene = sum(sig.all.gene), sig.up.gene = sum(sig.up.gene),
                                                              sig.down.gene = sum(sig.down.gene), 
                                                              alt_splice = sum(match))

gene.counts$isoform <- ifelse(gene.counts$threshold.all > 0, TRUE, FALSE)

gene.counts$gene <- ifelse(gene.counts$sig.all.gene > 0, TRUE, FALSE)

ggplot(results, aes(A = threshold.up, B = threshold.down)) +
  geom_venn()

# venn diagram at gene level

head(gene.counts)

ggplot(gene.counts, aes(A = isoform, B = gene)) +
  geom_venn()

# list of isoforms without corresponding significant genes

isoforms.only.up <- results[!is.na(results$padj.up),]
isoforms.only.up <- results[results$sig.up == TRUE & results$sig.up.gene == FALSE,]
isoforms.only.up.gene.counts <- isoforms.only.up %>% group_by(Gene.EnsID) %>% summarise(total.isoforms = n())
isoforms.only.up.gene.counts <- isoforms.only.up.gene.counts[order(isoforms.only.up.gene.counts$total.isoforms, decreasing = T),]
isoforms.only.up <- results[results$Gene.EnsID %in% isoforms.only.up.gene.counts$Gene.EnsID,]
isoforms.only.up <- isoforms.only.up[isoforms.only.up$sig.up == TRUE,]
isoforms.only.up <- isoforms.only.up[order(isoforms.only.up$padj.up),]
head(isoforms.only.up)
     

isoforms.only.down <- results[!is.na(results$padj.down),]
isoforms.only.down <- results[results$sig.down == TRUE & results$sig.down.gene == FALSE,]
isoforms.only.down.gene.counts <- isoforms.only.down %>% group_by(Gene.EnsID) %>% summarise(total.isoforms = n())
isoforms.only.down.gene.counts <- isoforms.only.down.gene.counts[order(isoforms.only.down.gene.counts$total.isoforms, decreasing = T),]
isoforms.only.down <- results[results$Gene.EnsID %in% isoforms.only.down.gene.counts$Gene.EnsID,]
isoforms.only.down <- isoforms.only.down[isoforms.only.down$sig.down == TRUE,]
isoforms.only.down <- isoforms.only.down[order(isoforms.only.down$padj.down),]
head(isoforms.only.down)

isoforms.only.all <- results[!is.na(results$padj.all),]
isoforms.only.all <- results[results$sig.all == TRUE & results$sig.all.gene == FALSE,]
isoforms.only.all.gene.counts <- isoforms.only.all %>% group_by(Gene.EnsID) %>% summarise(total.isoforms = n())
isoforms.only.all.gene.counts <- isoforms.only.all.gene.counts[order(isoforms.only.all.gene.counts$total.isoforms, decreasing = T),]
isoforms.only.all <- results[results$Gene.EnsID %in% isoforms.only.all.gene.counts$Gene.EnsID,]
isoforms.only.all <- isoforms.only.all[isoforms.only.all$sig.all == TRUE,]
isoforms.only.all <- isoforms.only.all[order(isoforms.only.all$padj.all),]
head(isoforms.only.all)

# add isoform only, no significant genes results to results data frame

results$isoform.only.all <- results$Transcript.EnsID %in% isoforms.only.all$Transcript.EnsID
results$isoform.only.up <- results$Transcript.EnsID %in% isoforms.only.up$Transcript.EnsID
results$isoform.only.down <- results$Transcript.EnsID %in% isoforms.only.down$Transcript.EnsID

ggplot(results, aes(A = isoform.only.all, B = isoform.only.up, C = isoform.only.down)) +
  geom_venn()


results[results$isoform.only.all == TRUE & results$isoform.only.up == TRUE & results$isoform.only.down,]

# evidence of alternative splicing between responder types at recurrence - 92 genes, corresponding to 

gene.counts.filter <- gene.counts[gene.counts$alt_splice > 1 & gene.counts$threshold.up > 0 & gene.counts$threshold.down > 0 &
                                    gene.counts$sig.up.gene == 0 & gene.counts$sig.down.gene == 0,]

isoforms.filter <- results[!is.na(results$Gene.EnsID),]
isoforms.filter <- isoforms.filter[(isoforms.filter$Gene.EnsID %in% gene.counts.filter$Gene.EnsID),]
isoforms.filter <- isoforms.filter[isoforms.filter$threshold.up == TRUE | isoforms.filter$threshold.down == TRUE,]

unique(sub("-.*","",isoforms.filter$GeneName))

isoforms.filter <- as.data.frame(isoforms.filter)

alt_splice <- isoforms.filter %>% dplyr::count(Gene.EnsID.Simp, sort = TRUE) %>% filter(n > 1)

isoforms.filter <- isoforms.filter[isoforms.filter$Gene.EnsID.Simp %in% alt_splice$Gene.EnsID.Simp,]
isoforms.filter <- isoforms.filter[isoforms.filter$threshold.up != isoforms.filter$threshold.down,]



                         
isoforms.filter <- isoforms.filter[order(isoforms.filter$padj.up),]
isoforms.filter[grep("EGFR",isoforms.filter$GeneName),]

results[grep("SNHG14",results$GeneName),]

test <- isoforms.filter %>% group_by(Gene.EnsID) %>% summarise(n())

table(test[,2])

isoforms.filter

isoforms.filter[isoforms.filter$jarid2.tss == TRUE,]
