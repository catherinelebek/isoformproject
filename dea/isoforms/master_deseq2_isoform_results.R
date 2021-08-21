# This script is primarily used to create volcano plots of the DEA results. This includes annotating DEIs by:
# Direction of dysregulation per responder-type
# Validation in GLASS dataset
# Relationship to gene-level DEA
# Relationship to JARID2 regulation

# The main results data frame used in this script is called "merge"


# import libraries ####

library(DESeq2) 
library(ggplot2)
library(ggrepel)
library(ggsci)
library(dplyr)
library(clipr)

# Fill in this before starting ####

# Is this script running for all patients, up-responders or down-responders?

patienttype <- "Up responders"
deseq2output <- "up-responders/"

# creating a results dataframe from DESeq2 results ######

load(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/",deseq2output,"deseq2.RData"))

# load full list of transcripts in order to pull through gene names

genelist <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

# extract just transcripts and gene names

genelist <- genelist[,c(1,10)]

# create results table from dds object

res <- results(dds, alpha = 0.05)

# put in a dataframe so can merge with the list of transcripts and gene names

resOrdered <- as.data.frame(res)

# run merge to pull through gene names into the DESeq2 results table

merge <- merge(resOrdered, genelist, by.x = "row.names", by.y = "EnsID", all.x = TRUE)

# order by p-value

merge <- merge[order(merge$pvalue),]

# move gene list to be second column

merge <- merge[,c(1,8,2:7)]

# get a summary of the results at the 90% confidence level

summary(res)

# check for isoforms with NAs for adjusted p-values

table(!is.na(merge$padj))
# merge <- merge[!is.na(merge$padj),]

# save results to csv

write.csv(merge, paste0("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/",deseq2output,"deseq2results.csv"))

# pull in mapping and jarid2 data #####

mapping <- read.delim("~/Documents/Semester3/Project/InputData/output.txt", header = F, sep = "\t") # gene-to-transcript IDs

colnames(mapping) <- c("Gene.EnsID","Transcript.EnsID")

mapping$Gene.EnsID.Simp <- sub("\\..*","",mapping$Gene.EnsID) # add a non-versioned gene Ens ID
mapping$Transcript.EnsID.Simp <- sub("\\..*","",mapping$Transcript.EnsID) # add a non-versioned transcript Ens ID

jarid2.gene <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",") # JARID2 gene IDs
jarid2.gene <- c(t(jarid2.gene))

jarid2.tss <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/jarid.tss.transcripts.csv",
                       header = T) # JARID2 TSS transcript IDs

jarid2.tss <- c(t(jarid2.tss))

jarid2.tss <- sub("\\..*","",jarid2.tss) # create a non-versioned list of transcript Ens IDs

# creating data frame for annotating transcripts ####

annotation <- mapping
annotation$jarid2.gene <- annotation$Gene.EnsID.Simp %in% jarid2.gene

table(annotation$jarid2.gene) # TRUE = number of transcripts in mapping data from JARID2 genes

annotation$jarid2.tss <- annotation$Transcript.EnsID.Simp %in% jarid2.tss

table(annotation$jarid2.tss) # TRUE = number of transcripts in mapping data with JARID2 TSSs

# merge results with annotations ####

head(annotation)
head(merge)

table(!merge$Row.names %in% annotation$Transcript.EnsID) # number of isoforms without gene annotations
merge[!merge$Row.names %in% annotation$Transcript.EnsID,] # list of isoforms without gene annotations, ordered by adjusted p-value

merge <- left_join(merge, annotation, by = c("Row.names" = "Transcript.EnsID")) # merge datasets

head(merge)

sum(is.na(merge$Gene.EnsID)) # rows that are in the results table but not in the annotation
table(is.na(merge$Gene.EnsID) & merge$padj < 0.05) # check if any are significant

# create columns for annotating a volcano plot ####

merge$threshold <- merge$padj < 0.05 & abs(merge$log2FoldChange) > 1

merge$jarid.status <- ifelse(is.na(merge$jarid2.gene), "Not in mapping",
                             ifelse(merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE, "JARID2 TSS",
                             ifelse(merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE, "JARID2 Gene",
                             "Not JARID2-regulated")))
table(merge$jarid.status)

merge$Legend <- ifelse(merge$threshold == FALSE, "Below threshold", merge$jarid.status)

# format columns as factors to aid in creating volcano plots

merge$jarid.status <- as.factor(merge$jarid.status)
merge$Legend <- as.factor(merge$Legend)
merge$threshold <- as.factor(merge$threshold)

# rerrange columns ####

merge <- merge[,c(1,11,9,10,2,3:8,12:16)]
head(merge)

# import data on isoforms only signficant in the combined DEA vs. split by responder-type ####

justall <- read.table("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/justall.txt", header = F)
justall <- c(t(justall))

merge$justall <- merge$Row.names %in% justall

table(merge$justall)

# save results to csv ####

write.csv(merge, paste0("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/",deseq2output,"deseq2results_merge.csv"))


# create volcano plot of overall results, colouring up- and down-regulated isoforms above the threshold differently ####

merge$top <- ifelse((-log10(merge$padj) > 4 & merge$log2FoldChange > 1) | 
                      (-log10(merge$padj) > 9 & merge$log2FoldChange < -1) ,1,0) # column to use for selecting data points to label

merge$direction <- ifelse(merge$threshold == TRUE & merge$log2FoldChange > 1, "Up-regulated",
                          ifelse(merge$threshold == TRUE & merge$log2FoldChange < -1, "Down-regulated",
                                 "Below Threshold")) # column to use for colouring datapoints

p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = direction, shape = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, top == 1),
                   aes(label = GeneName, fill = direction),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  scale_color_manual(name = " ", values = c("Up-regulated" = "red",
                                "Down-regulated" = "blue",
                                "Below threshold" = "grey")) +
  scale_fill_manual(values = c("Up-regulated" = "red",
                               "Down-regulated" = "blue",
                               "Below threshold" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "direction") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(1,0,0.5,0.5), "cm"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        axis.text.x = element_text(size = "12"),
        axis.text.y = element_text(size = "12"),
        legend.title = element_text(face = "bold", size = "12"),
        legend.text = element_text(size = "12")) +
  ggpubr::rremove("grid")


p

# create volcano plot highlighting genes by JARID2 relationship ####

merge$top <- ifelse((-log10(merge$padj) > 17 & merge$log2FoldChange > 1) | 
                      (-log10(merge$padj) > 9 & merge$log2FoldChange < -1) ,1,0)


p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Legend, shape = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, top == 1),
                   aes(label = GeneName, fill = Legend),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  scale_color_manual(name = " ", values = c("JARID2 Gene" = "#0609F1", 
                                "JARID2 TSS" = "#F11406",
                                "Not JARID2-regulated" = "#208C17",
                                "Below threshold" = "grey")) +
  scale_fill_manual(values = c("JARID2 Gene" = "#0609F1", 
                                "JARID2 TSS" = "#F11406",
                                "Not JARID2-regulated" = "#208C17",
                                "Below threshold" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Legend") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(1,0,0.5,0.5), "cm"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        axis.text.x = element_text(size = "12"),
        axis.text.y = element_text(size = "12"),
        legend.title = element_text(face = "bold", size = "12"),
        legend.text = element_text(size = "12")) +
  ggpubr::rremove("grid")


p

# pull in gene level significance data ####

gene.dea.res <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/genes/",deseq2output,"deseq2results.csv"),
                           header = T, sep = ",")

gene.dea.res <- gene.dea.res[,c(2,5,9)] # only keep gene Ens IDs, L2FC & padj columns

gene.dea.res.sig <- gene.dea.res[gene.dea.res$padj < 0.05,] # filter for significant genes

# annotate DEIs based on whether they have a corresponding DEG

merge$gene.status <- merge$Gene.EnsID %in% gene.dea.res.sig$Row.names
merge$gene.status <- ifelse(merge$threshold == TRUE & merge$gene.status == FALSE, "No DEG", 
                            ifelse(merge$threshold == TRUE & merge$gene.status == TRUE, "DEG",
                                   "Below Threshold"))

# combine merge data frame and gene-level results

merge <- merge(merge, gene.dea.res, by.x = "Gene.EnsID", by.y = "Row.names", all.x = TRUE)

# annotate DEIs with associated DEGs based on whether the direction of dysregulation is the same

table(merge$gene.status)

merge$gene.isoform.direction <- ifelse(merge$gene.status == "DEG" &
                                        merge$log2FoldChange.x*merge$log2FoldChange.y < 0,
                                       "DEG - different direction", 
                                       ifelse(merge$gene.status == "DEG", "DEG - same direction",
                                              merge$gene.status))

# update column for labelling data points in volcano plot

merge$top <- ifelse((merge$gene.isoform.direction == "DEG - different direction" &
                       -log10(merge$padj.x) > 1) |
                    (merge$gene.isoform.direction == "No DEG" & 
                      -log10(merge$padj.x) > 6) , 1, 0)

p <- ggplot(merge, aes(x = log2FoldChange.x, y = -log10(padj.x))) +
  geom_point(aes(colour = gene.isoform.direction, shape = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, top == 1),
                   aes(label = GeneName, fill = gene.isoform.direction),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  scale_color_manual(values = c("DEG - different direction" = "darkgreen",
                                "No DEG" = "#F11406", 
                                "DEG - same direction" = "black",
                                "Below Threshold" = "grey")) + 
  scale_fill_manual(values =  c("DEG - different direction" = "darkgreen",
                                "No DEG" = "#F11406", 
                                "DEG - same direction" = "black",
                                "Below Threshold" = "grey")) + 
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Gene Status") +
  theme_bw() +
  theme(plot.margin = unit(c(1,0,0.5,0.5), "cm"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        axis.text.x = element_text(size = "12"),
        axis.text.y = element_text(size = "12"),
        legend.title = element_text(face = "bold", size = "12"),
        legend.text = element_text(size = "12")) +
  ggpubr::rremove("grid")

p

table(merge$gene.isoform.direction) # particularly interested in DEIs with DEGs in different direction - 34


# compare to glass results ######

glass.up <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/glass/",deseq2output,"glassfilter/deseq2results.csv"),
                       sep = ",", header = T)


head(merge)
head(glass.up)

# combine merge and glass datasets using non-versioned isoform Ens ID (as version not available for GLASS data)
# this step is a little bit slow (takes a couple of minutes)

merge <- merge(merge, glass.up[,c("X","padj","log2FoldChange","pvalue")], by.x = "Transcript.EnsID.Simp", by.y = "X", all.x = TRUE)

head(merge)

merge$glass.up <- ifelse(merge$padj < 0.05 & merge$log2FoldChange.x*merge$log2FoldChange > 0, "Both significant, same direction dysregulation", 
                         ifelse(merge$padj < 0.05, "Both significant, different direction dysregulation", "Not significant in GLASS"))

merge$glass.up <- ifelse(is.na(merge$glass.up), "Not analysed in GLASS", merge$glass.up)

merge$glass.up <- ifelse(merge$threshold == FALSE, "Below Threshold", merge$glass.up)

merge$glass.up <- as.factor(merge$glass.up)


table(merge$glass.up)

merge$top <- ifelse(merge$glass.up == "Both significant, same direction dysregulation" |
                    -log10(merge$padj.x) > 10, 1, 0)

p <- ggplot(merge, aes(x = log2FoldChange.x, y = -log10(padj.x))) +
  geom_point(aes(colour = glass.up, shape = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, top == 1),
                   aes(label = GeneName, fill = glass.up),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  scale_color_manual(values = c("Both significant, same direction dysregulation" = "darkgreen", 
                                "Not analysed in GLASS" = "black",
                                "Not significant in GLASS" = "red",
                                "Below Threshold" = "grey")) +
  scale_fill_manual(values = c("Both significant, same direction dysregulation" = "darkgreen", 
                               "Not analysed in GLASS" = "black",
                               "Not significant in GLASS" = "red",
                               "Below Threshold" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Comparison with GLASS data") +
  theme_bw() +
  theme(plot.margin = unit(c(1,0,0.5,0.5), "cm"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        axis.text.x = element_text(size = "12"),
        axis.text.y = element_text(size = "12"),
        legend.title = element_text(face = "bold", size = "12"),
        legend.text = element_text(size = "12")) +
  ggpubr::rremove("grid")

p


# comparing with glass data #####

test1 <- merge[!is.na(merge$padj.x),]
test1 <- test1[test1$padj.x < 0.05,]
test1 <- test1[!(is.na(test1$padj)),]
nrow(test1[test1$padj < 0.05,])


test2 <- merge[!is.na(merge$padj.x),]
test2 <- test2[test2$padj.x < 0.05 & abs(test2$log2FoldChange.x) > 1,]
test2 <- test2[!(is.na(test2$padj.y)),]
nrow(test2[test2$padj.y < 0.05 & abs(test2$log2FoldChange.y) > 1,])

table(merge$padj.x < 0.05)
table(merge$padj.x > 0.05 & merge$padj.y < 0.05)
table(merge$padj.x < 0.05 & abs(merge$log2FoldChange.x) > 1 & merge$padj.y < 0.05)
table(!is.na(merge$padj.x) & !is.na(merge$padj.y))


plot(test$padj.x, test$padj.y)
plot(test$padj.x, test$padj.y, xlim = c(0,1), ylim = c(0,1))
hist(test$padj.y, breaks = c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
cor(test1$padj.x, test1$padj)

cor(test$pvalue.x, test$pvalue.y)

cont.table <- data.frame("V1" = c(18,1364), "V2" = c(238,45476))
cont.table

chisq.test(cont.table)

# volcano plot of isoforms only found in all ####
# only run this section if all patient data being used

merge$justall <- ifelse(merge$threshold == FALSE, "Below Threshold", 
                        ifelse(merge$justall == FALSE, "Significant",
                               "Significant in all patient analysis only"))

table(merge$justall)
merge$top <- ifelse(merge$justall == "Significant in all patient analysis only", 1, 0)

p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = justall, shape = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, top == 1),
                   aes(label = GeneName, fill = justall),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  scale_color_manual(name = " ", values = c("Significant in all patient analysis only" = "red",
                                            "Significant" = "#208C17",
                                            "Below Threshold" = "grey")) +
  scale_fill_manual(values = c("Significant in all patient analysis only" = "red",
                                "Significant" = "#208C17",
                                "Below Threshold" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "direction") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = unit(c(1,0,0.5,0.5), "cm"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        axis.text.x = element_text(size = "12"),
        axis.text.y = element_text(size = "12"),
        legend.title = element_text(face = "bold", size = "12"),
        legend.text = element_text(size = "12")) +
  ggpubr::rremove("grid")


p


        
