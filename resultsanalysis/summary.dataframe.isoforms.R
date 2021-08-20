# This script is used to summarise all the results of the paired DEA for the in-house data
# The results are summarised at the isoform level in a dataframe called "results"
# This dataframe has 218,712 rows, corresponding to all the isoforms quantified

# import libraries ####

library(tidyverse) # for manipulating data frames
library(VennDiagram) # for plotting venn diagrams
library(ggvenn) # for formatting venn diagrams
library(clipr) # for copying data objects out of R Studio directly


# full list of 218,712 isoforms ####

isoforms <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")
 
isoforms <- isoforms[,c(1,10)] # retaining just the isoform EnsID and gene name columns

# add in mapping data to annotate transcripts with their corresponding gene ID and create results dataframe ####

mapping <- read.delim("~/Documents/Semester3/Project/InputData/output.txt", header = F, sep = "\t") # gene-to-transcript IDs

colnames(mapping) <- c("Gene.EnsID","Transcript.EnsID")

mapping$Gene.EnsID.Simp <- sub("\\..*","",mapping$Gene.EnsID) # adding non-versioned gene Ens ID to mapping data
mapping$Transcript.EnsID.Simp <- sub("\\..*","",mapping$Transcript.EnsID) # adding non-versioned isoform Ens ID to mapping data

table(isoforms$EnsID %in% mapping$Transcript.EnsID) # check how many of the 218,712 isoforms appear in the mapping data

# creating results data frame

results <- left_join(isoforms, mapping, by = c("EnsID" = "Transcript.EnsID"))

# import isoform paired DEA reuslts ####

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

# make sure all the isoforms are unique

nrow(results)
length(unique(results$Transcript.EnsID))

nrow(isoforms.all)
length(unique(isoforms.all$Row.names))

nrow(isoforms.up)
length(unique(isoforms.up$Row.names))

nrow(isoforms.down)
length(unique(isoforms.down$Row.names))

# adding which isoforms were tested in each DEA ####

results$tested.all <- isoforms$EnsID %in% isoforms.all$Row.names
results$tested.up <- isoforms$EnsID %in% isoforms.up$Row.names
results$tested.down <- isoforms$EnsID %in% isoforms.down$Row.names

head(results)

colnames(results)[1] <- "Transcript.EnsID" # rename first column to specify isoform (not gene) Ens ID

# importJARID2 data ####

jarid2.gene <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                          header = F, sep = ",") # JARID2 gene IDs

jarid2.gene <- c(t(jarid2.gene))

jarid2.tss <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/jarid.tss.transcripts.csv",
                       header = T) # JARID2 TSS transcript IDs

jarid2.tss <- c(t(jarid2.tss))
jarid2.tss <- sub("\\..*","",jarid2.tss)

# add JARID2 data to results data frame ####

results$jarid2.gene <- results$Gene.EnsID.Simp %in% jarid2.gene # joining on non-versioned gene Ens ID

results$jarid2.tss <- results$Transcript.EnsID.Simp %in% jarid2.tss # joining on non-versioned isoform Ens ID

head(results)

# add isoform DEA results to the main results data frame ####

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

# import the gene-level DEA data ####

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

# add gene DEA results to the main results data frame ####

results <- merge(results, genes.all[,c("EnsID","log2FoldChange","padj")], by.x = "Gene.EnsID", by.y = "EnsID", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.all.genes","padj.all.genes")


results <- merge(results, genes.up[,c("EnsID","log2FoldChange","padj")], by.x = "Gene.EnsID", by.y = "EnsID", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.up.genes","padj.up.genes")


results <- merge(results, genes.down[,c("EnsID","log2FoldChange","padj")], by.x = "Gene.EnsID", by.y = "EnsID", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.down.genes","padj.down.genes")

head(results)

# rearrange columns ####

results <- results[,c(2,5,3,1,4,9,10,6,7,8,11:22)]
head(results)

# add columns to indicate which isoforms DE in each isoform-level DEA ####

results$sig.all <- ifelse(is.na(results$padj.all), FALSE, ifelse(results$padj.all < 0.05, TRUE, FALSE))
results$sig.up <- ifelse(is.na(results$padj.up), FALSE, ifelse(results$padj.up < 0.05, TRUE, FALSE))
results$sig.down <- ifelse(is.na(results$padj.down), FALSE, ifelse(results$padj.down < 0.05, TRUE, FALSE))

# column to indicate if isoform tested in both up and down responder DEA

results$tested.up.down <- ifelse(results$tested.up == TRUE & results$tested.down == TRUE, TRUE, FALSE)

results$threshold.all <- results$sig.all == TRUE & abs(results$LFC.all) > 1
results$threshold.up <- results$sig.up == TRUE & abs(results$LFC.up) > 1
results$threshold.down <- results$sig.down == TRUE & abs(results$LFC.down) > 1

# check results as expected

table(results$sig.all)
table(results$sig.up)
table(results$sig.down)

table(results$threshold.all)
table(results$threshold.up)
table(results$threshold.down)


# add columns to indicate which genes DE in each gene-level DEA ####

results$sig.all.gene <- ifelse(is.na(results$padj.all.genes), FALSE, ifelse(results$padj.all.genes < 0.05, TRUE, FALSE))
results$sig.up.gene <- ifelse(is.na(results$padj.up.genes), FALSE, ifelse(results$padj.up.genes < 0.05, TRUE, FALSE))
results$sig.down.gene <- ifelse(is.na(results$padj.down.genes), FALSE, ifelse(results$padj.down.genes < 0.05, TRUE, FALSE))

# check results as expected

table(results$sig.all.gene)
table(results$sig.up.gene)
table(results$sig.down.gene)

# Venn Diagrams ####

# isoforms tested - not included in final report

ggplot(results, aes(A = tested.up, B = tested.down)) +
  geom_venn(set_names = c("Up-responders","Down-responders")) +
  theme_bw() +
  ggpubr::rremove("grid") +
  ggpubr::rremove("axis")

# jarid2 genes vs. jarid2 genes & tss - not included in final report

ggplot(results, aes(A = jarid2.gene, B = jarid2.tss)) +
  geom_venn()

# isoforms tested in both analyses and/or significant across up and/or down-responders - included in final report

ggplot(results, aes(A = sig.up, B = sig.down, C = tested.up.down)) +
  geom_venn(set_names = c("Up",
                          "Down",
                          "Both")) +
  theme_bw() +
  ggpubr::rremove("grid") +
  ggpubr::rremove("axis")

# isoforms in up and down-responders that are above the threshold - not included in final report

ggplot(results, aes(A = threshold.up, B = threshold.down)) +
  geom_venn()

# bar plot showing proportion of DEIs found in both responder types that are dysregulated in the same direction ####

sig.results.both <- results[results$sig.up == TRUE & results$sig.down == TRUE,] # signficant resuls in both
sig.results.both$direction <- ifelse(sig.results.both$LFC.up*sig.results.both$LFC.down > 0, "Same", "Different")

p <- ggplot(data = sig.results.both, aes(direction)) +
  geom_bar() +
  theme_bw() +
  scale_x_discrete(name = "Direction", limits = c("Same","Different")) +
  theme(axis.text.x = element_text(size = "12"), 
        axis.text.y = element_text(size = "12"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  ggpubr::rremove("grid")

p

clipr(sig.results.both[sig.results.both$direction == "Same",]) # copy out table of 15 isoforms that are dysregulated in the same direction

# summarise at the gene level ####

gene.counts <- results %>% group_by(Gene.EnsID) %>% summarise(total.isoforms = n(), threshold.all = sum(threshold.all), 
                                                              threshold.up = sum(threshold.up), threshold.down = sum(threshold.down),
                                                              sig.all = sum(sig.all), sig.up = sum(sig.up), sig.down = sum(sig.down),
                                                              sig.all.gene = sum(sig.all.gene), sig.up.gene = sum(sig.up.gene),
                                                              sig.down.gene = sum(sig.down.gene), sig.up.up = sum(LFC.up > 0, na.rm = TRUE & sig.up),
                                                              sig.up.down = sum(sig.up & LFC.up < 0, na.rm = TRUE), sig.down.down = sum(sig.down & LFC.down < 0, na.rm = TRUE),
                                                              sig.down.up = sum(sig.down & LFC.down > 0, na.rm = TRUE), sig.up.up.gene = sum(sig.up.gene & LFC.up.genes > 0),
                                                              sig.up.down.gene = sum(sig.up.gene & LFC.up.genes < 0), sig.down.down.gene = sum(sig.down.gene & LFC.down.genes < 0),
                                                              sig.down.up.gene = sum(sig.down.gene & LFC.down.genes > 0))


gene.counts <- as.data.frame(gene.counts)
write.csv(gene.counts, "~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/genecounts.csv")

# isoforms that are significant/more significant in all patients compared to splitting out responders separately ####

all <- results[results$tested.up == TRUE & results$padj.all < results$padj.up & results$padj.all < results$padj.down,]
all <- all[!is.na(all$padj.all),]
all <- all[all$padj.all < 0.05,]
head(all)

all <- all[order(all$padj.all),]
justall <- all[all$sig.up == FALSE & all$sig.down == FALSE,1]

write.table(justall, "~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/justall.txt", row.names = FALSE, col.names = F)

# pull in GLASS data

glass.up <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/glass/up-responders/glassfilter/deseq2results.csv",
                                    header = T, sep = ",")

colnames(glass.up)[1] <- "EnsID"

glass.down <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/glass/down-responders/glassfilter/deseq2results.csv",
                       header = T, sep = ",")

colnames(glass.down)[1] <- "EnsID"

results <- merge(results, glass.up[,c("EnsID","log2FoldChange","padj")], by.x = "Transcript.EnsID.Simp", by.y = "EnsID", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.glass.up","padj.glass.up")

results <- merge(results, glass.down[,c("EnsID","log2FoldChange","padj")], by.x = "Transcript.EnsID.Simp", by.y = "EnsID", all.x = TRUE)
n <- ncol(results)
colnames(results)[(n-1):n] <- c("LFC.glass.down","padj.glass.down")

head(results)



knitr::write_bib(, "~/Documents/Semester3/Project/Report/citations/R.bib")