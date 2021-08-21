# import packages #### 

library(tidyverse)
library(clipr)

# import all data ####

isoforms.1 <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv", header = T, sep = ",")
isoforms.2 <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv", header = T, sep = ",")
isoforms.3 <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results.csv", header = T, sep = ",")
summary <- read.delim("~/Documents/Semester3/Project/Results/resultsanalysis/summarydf.csv", header = T, sep = ",")

# Top 20 isoforms for all patients as a starting point ####

isoforms <- isoforms.1 %>% select(Row.names , GeneName, log2FoldChange, padj)
isoforms <- dplyr::rename(isoforms, EnsID = Row.names)

merge <- isoforms.2 %>% select(Row.names, log2FoldChange, padj)
isoforms <- left_join(isoforms, merge, by = c(EnsID = "Row.names"), suffix = c(".all",".up"))

merge <- isoforms.3 %>% select(Row.names, log2FoldChange, padj)
isoforms <- left_join(isoforms, merge, by = c(EnsID = "Row.names"))

isoforms <- rename(isoforms, c(log2FoldChange.down = log2FoldChange, padj.down = padj))
isoforms$GeneNameSimp <- sub("-.*","",isoforms$GeneName)


merge <- summary %>% select(GeneName, jarid2, pca.gene.rank, padj.genes.all, LFC.genes.all, 
                            padj.genes.up, LFC.genes.up, padj.genes.down, LFC.genes.down, isoforms.full.counts, isoforms.all.counts.sig, isoforms.up.counts.sig,
                            isoforms.down.counts.sig)

isoforms <- left_join(isoforms, merge, by = c("GeneNameSimp" = "GeneName"))

top <- head(isoforms,20)

write_clip(top)


top


# Top 20 isoforms for up-responders as a starting point ####

isoforms <- isoforms.2 %>% select(Row.names , GeneName, log2FoldChange, padj)
isoforms <- dplyr::rename(isoforms, EnsID = Row.names)

merge <- isoforms.1 %>% select(Row.names, log2FoldChange, padj)
isoforms <- left_join(isoforms, merge, by = c(EnsID = "Row.names"), suffix = c(".all",".up"))

merge <- isoforms.3 %>% select(Row.names, log2FoldChange, padj)
isoforms <- left_join(isoforms, merge, by = c(EnsID = "Row.names"))

isoforms <- rename(isoforms, c(log2FoldChange.down = log2FoldChange, padj.down = padj))
isoforms$GeneNameSimp <- sub("-.*","",isoforms$GeneName)


merge <- summary %>% select(GeneName, jarid2, pca.gene.rank, padj.genes.all, LFC.genes.all, 
                            padj.genes.up, LFC.genes.up, padj.genes.down, LFC.genes.down, isoforms.full.counts, isoforms.all.counts.sig, isoforms.up.counts.sig,
                            isoforms.down.counts.sig)

isoforms <- left_join(isoforms, merge, by = c("GeneNameSimp" = "GeneName"))

top <- head(isoforms,20)

top

write_clip(top)

head(summary)
# Top 20 isoforms for down-responders as a starting point ####

isoforms <- isoforms.3 %>% select(Row.names , GeneName, log2FoldChange, padj)
isoforms <- dplyr::rename(isoforms, EnsID = Row.names)

merge <- isoforms.1 %>% select(Row.names, log2FoldChange, padj)
isoforms <- left_join(isoforms, merge, by = c(EnsID = "Row.names"), suffix = c(".all",".up"))

merge <- isoforms.2 %>% select(Row.names, log2FoldChange, padj)
isoforms <- left_join(isoforms, merge, by = c(EnsID = "Row.names"))

isoforms <- rename(isoforms, c(log2FoldChange.down = log2FoldChange, padj.down = padj))
isoforms$GeneNameSimp <- sub("-.*","",isoforms$GeneName)


merge <- summary %>% select(GeneName, jarid2, pca.gene.rank, padj.genes.all, LFC.genes.all, 
                            padj.genes.up, LFC.genes.up, padj.genes.down, LFC.genes.down, isoforms.full.counts, isoforms.all.counts.sig, isoforms.up.counts.sig,
                            isoforms.down.counts.sig)

isoforms <- left_join(isoforms, merge, by = c("GeneNameSimp" = "GeneName"))

top <- head(isoforms,20)

top

write_clip(top)

# other ####

head(isoforms.2,10)

isoforms.2[sub("-.*","",isoforms.2$GeneName) == "JARID2",]
isoforms.3[sub("-.*","",isoforms.3$GeneName) == "JARID2",]


head(summary)
head(isoforms.2)


upsig <- isoforms.2[isoforms.2$padj < 0.05,3]

upsig <- sub("-.*","",upsig)
upsig

write_clip(upsig)

summary[summary$GeneName == "JARID2",]
