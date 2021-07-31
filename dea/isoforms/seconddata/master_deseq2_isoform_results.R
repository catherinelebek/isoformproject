# import libraries ####

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(ggsci)
library(clipr)

# Fill in this before starting ####

# Is this script running for all patients, up-responders or down-responders? #### All Patients

patienttype <- ""
# deseq2output <- "" # All patients
deseq2output <- "down-responders/"
# deseq2output <- "down-responders"


# creating a results dataframe ######

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

# remove NAs

merge <- merge[!is.na(merge$padj),]

# save results to csv

write.csv(merge, paste0("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/",deseq2output,"deseq2results.csv"))

# pull in mapping and jarid2 data #####

mapping <- read.delim("~/Documents/Semester3/Project/InputData/output.txt", header = F, sep = "\t") # gene-to-transcript IDs

colnames(mapping) <- c("Gene.EnsID","Transcript.EnsID")

mapping$Gene.EnsID.Simp <- sub("\\..*","",mapping$Gene.EnsID)
mapping$Transcript.EnsID.Simp <- sub("\\..*","",mapping$Transcript.EnsID)

jarid2.gene <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",") # JARID2 gene IDs
jarid2.gene <- c(t(jarid2.gene))

jarid2.tss <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/jarid.tss.transcripts.csv",
                       header = T) # JARID2 TSS transcript IDs

jarid2.tss <- c(t(jarid2.tss))
jarid2.tss <- sub("\\..*","",jarid2.tss)


# creating data frame for annotating transcripts ####

annotation <- mapping
annotation$jarid2.gene <- annotation$Gene.EnsID.Simp %in% jarid2.gene

table(annotation$jarid2.gene)

annotation$jarid2.tss <- annotation$Transcript.EnsID.Simp %in% jarid2.tss

table(annotation$jarid2.tss)

# merge results with annotations ####

head(annotation)
head(merge)

table(!merge$Row.names %in% annotation$Transcript.EnsID) # number of isoforms without gene annotations
merge[!merge$Row.names %in% annotation$Transcript.EnsID,] # list of isoforms without gene annotations, ordered by adjusted p-value

merge <- left_join(merge, annotation, by = c("Row.names" = "Transcript.EnsID"))

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

merge$top <- ifelse((-log10(merge$padj) > 15 & merge$log2FoldChange > 1) | 
                     (merge$padj < 0.05 & merge$log2FoldChange > 10) |
                     (-log10(merge$padj) > 10 & merge$log2FoldChange < -1) |
                      (merge$log2FoldChange < -5),1,0)

merge$Legend <- ifelse(merge$threshold == FALSE, "Below threshold", merge$jarid.status)

# save results to csv

write.csv(merge, paste0("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/",deseq2output,"deseq2results_merge.csv"))

merge$jarid.status <- as.factor(merge$jarid.status)
merge$Legend <- as.factor(merge$Legend)
merge$threshold <- as.factor(merge$threshold)

merge$order <- ifelse(merge$jarid.status == "JARID2 TSS", 1, 
                      ifelse(merge$jarid.status == "JARID2 Gene", 2,
                      3))

merge <- merge[order(merge$order, decreasing = T),]
head(merge)

table(merge$order)

# create volcano plot ####

p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Legend, shape = threshold)) +
  ggtitle(paste0("Differential Isoform Expression - ", patienttype)) +
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
  scale_color_manual(values = c("JARID2 Gene" = "#0609F1", 
                                "JARID2 TSS" = "#F11406",
                                "Not JARID2-regulated" = "#72C97E",
                                "Below threshold" = "grey")) +
  scale_fill_manual(values = c("JARID2 Gene" = "#0609F1", 
                                "JARID2 TSS" = "#F11406",
                                "Not JARID2-regulated" = "#72C97E",
                                "Below threshold" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Legend") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p

# plot just JARID2 genes ####

merge$top <- ifelse(((-log10(merge$padj)) > 4) | 
                    (abs(merge$log2FoldChange) > 10 & (-log10(merge$padj) > 1.4)),
                     1,0)

p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = jarid.status, shape = threshold)) +
  ggtitle(paste0("Differential Isoform Expression - ", patienttype)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, jarid2.gene == TRUE & threshold == TRUE),
                   aes(label = GeneName, fill = jarid.status),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0, "lines"),
                   segment.colour = "black") +
  scale_color_manual(values = c("JARID2 Gene" = "#0609F1", 
                                "JARID2 TSS" = "#F11406",
                                "Not JARID2-regulated" = "grey")) +
  scale_fill_manual(values = c("JARID2 Gene" = "#0609F1", 
                                "JARID2 TSS" = "#F11406",
                                "Not JARID2-regulated" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "JARID2 Status") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p

# plot just JARID2 TSS transcripts ####

merge$top <- ifelse(((-log10(merge$padj)) > 4) | 
                      (abs(merge$log2FoldChange) > 10 & (-log10(merge$padj) > 1.4)),
                    1,0)

p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = jarid.status, shape = threshold)) +
  ggtitle(paste0("Differential Isoform Expression - ", patienttype)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, jarid2.tss == TRUE & threshold == TRUE),
                   aes(label = GeneName, fill = jarid.status),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0, "lines"),
                   segment.colour = "black") +
  scale_color_manual(values = c("JARID2 Gene" = "#0609F1", 
                                "JARID2 TSS" = "#F11406",
                                "Not JARID2-regulated" = "grey")) +
  scale_fill_manual(values = c("JARID2 Gene" = "#0609F1", 
                               "JARID2 TSS" = "#F11406",
                               "Not JARID2-regulated" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "JARID2 Status") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p

# pull in gene level significance data ####

gene.dea.res <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/genes/",deseq2output,"/deseq2results.csv"),
                           header = T, sep = ",")

gene.dea.res <- gene.dea.res[,c(2,5,9)]

gene.dea.res.sig <- gene.dea.res[gene.dea.res$padj < 0.05,]

merge$gene.status <- merge$Gene.EnsID %in% gene.dea.res.sig$Row.names

p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = gene.status, shape = threshold)) +
  ggtitle(paste0("Differential Isoform Expression - ", patienttype)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, gene.status == FALSE & threshold == TRUE),
                   aes(label = GeneName, fill = gene.status),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.05, "lines"),
                   segment.colour = "black") +
  scale_color_manual(values = c("TRUE" = "grey50", 
                                "FALSE" = "#F11406")) + 
  scale_fill_manual(values = c("TRUE" = "#0609F1", 
                               "FALSE" = "#F11406")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Gene Status") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p


# pull in up and down responder data - is the dysregulation in the same or opposite direction?

up.dea.res <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv"),
                           header = T, sep = ",")
down.dea.res <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results.csv"),
                           header = T, sep = ",")

head(up.dea.res)
head(down.dea.res)

up.dea.res <- up.dea.res[,c(2,3,5,9)]
down.dea.res <- down.dea.res[,c(2,3,5,9)]

up.dea.res.sig <- up.dea.res[up.dea.res$padj < 0.05,]
down.dea.res.sig <- down.dea.res[down.dea.res$padj < 0.05,]

merge <- left_join(merge, up.dea.res.sig[,c("Row.names","log2FoldChange")], 
                   by = c("Row.names" = "Row.names"), suffix = c("",".up"))


merge <- left_join(merge, down.dea.res.sig[,c("Row.names","log2FoldChange")], 
                   by = c("Row.names" = "Row.names"), suffix = c("",".down"))


merge$direction <- ifelse(is.na(merge$log2FoldChange.up) | is.na(merge$log2FoldChange.down),
                          "Unclear",
                           ifelse((merge$log2FoldChange.up > 0 & merge$log2FoldChange.down > 0) |
                                   (merge$log2FoldChange.up < 0 & merge$log2FoldChange.down < 0),
                          "Same","Opposite"))
table(merge$direction)

merge$direction <- as.factor(merge$direction)

merge$top <- ifelse(merge$direction == "Same", 1, 0)

merge <- merge[order(merge$top, decreasing = T),]

p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = direction, shape = threshold)) +
  ggtitle(paste0("Differential Isoform Expression - ", patienttype)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, direction == "Same" & threshold == TRUE),
                   aes(label = GeneName, fill = direction),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  scale_color_manual(values = c("Same" = "#F11406", 
                                "Opposite" = "blue",
                                "Unclear" = "grey")) + 
  scale_fill_manual(values = c("Same" = "#F11406", 
                               "Opposite" = "#blue",
                               "Unclear" = "yellow")) + 
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Direction") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p

# identifying transcripts for which a different transcript of the same gene was dysregulated in the other responder type ####

head(merge)
sum(!is.na(merge$log2FoldChange.down))

# get the genes which have transcripts significant in one responder type but not the other
potentialgenes <- sub("-.*","",merge[!is.na(merge$log2FoldChange.down) & is.na(merge$log2FoldChange.up), 2])

# make this a unique list
potentialgenes <- unique(potentialgenes)

# flag if a transcript in the main data frame is in this list of genes
merge$potentialgenes <- sub("-.*","",merge$GeneName) %in% potentialgenes

# if the transcript is in this list of genes then ask - is it significant in only one type of responder?
merge$potentialtranscripts <- (merge$potentialgenes == TRUE & !is.na(merge$log2FoldChange.up) & is.na(merge$log2FoldChange.down)) |
                              (merge$potentialgenes == TRUE & is.na(merge$log2FoldChange.up) & !is.na(merge$log2FoldChange.down))

head(merge[merge$potentialtranscripts == TRUE & !is.na(merge$log2FoldChange.up),])


merge$top <- ifelse((((-log10(merge$padj)) > 5) | 
                      (abs(merge$log2FoldChange) > 10 & (-log10(merge$padj) > 1.4)) &
                       merge$potentialtranscripts == TRUE), 1, 0)


table(merge$top)
# if the genes.alt.expr column value is TRUE, this means that the gene relating to the transcript has isoform switching differences between up- and down-responders

p <- ggplot(merge, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = potentialtranscripts, shape = threshold)) +
  ggtitle(paste0("Differential Isoform Expression - ", patienttype)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(merge, top == 1),
                   aes(label = GeneName, fill = potentialtranscripts),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  scale_color_manual(values = c("TRUE" = "#F11406", 
                                "FALSE" = "grey")) +
  scale_fill_manual(values = c("TRUE" = "#F11406", 
                               "FALSE" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Alternative Splicing") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p



#######


# compare to glass results ######

head(merge)

glass.up <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/glass/up-responders/deseq2results.csv",
                       sep = ",", header = T)

merge$EnsIDSimp <- sub("\\..*","",merge$Row.names)
head(merge)

merge <- merge(merge, glass.up[,c("X","padj","log2FoldChange")], by.x = "EnsIDSimp", by.y = "X", all.x = TRUE)

merge$glass.up <- ifelse(merge$padj.y < 0.05 & merge$log2FoldChange.x*merge$log2FoldChange.y > 0, "Both Significant, Same Direction Dysregulation", 
                         ifelse(merge$padj.y < 0.05, "Both Significant, Different Direction Dysregulation", "Different Results"))

merge$glass.up <- ifelse(is.na(merge$glass.up), "Not analysed in GLASS", merge$glass.up)

merge$glass.up <- ifelse(merge$padj.x > 0.05, "Below Threshold", merge$glass.up)

p <- ggplot(merge, aes(x = log2FoldChange.x, y = -log10(padj.x))) +
  geom_point(aes(colour = glass.up, shape = threshold)) +
  ggtitle(paste0("Differential Isoform Expression - ", patienttype)) +
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
  scale_color_manual(values = c("Both Significant, Same Direction Dysregulation" = "darkgreen", 
                                "Both Significant, Different Direction Dysregulation" = "orange",
                                "Not analysed in GLASS" = "black",
                                "Different Results" = "red",
                                "Below Threshold" = "grey")) +
  scale_fill_manual(values = c("Both Significant, Same Direction Dysregulation" = "darkgreen", 
                               "Both Significant, Different Direction Dysregulation" = "orange",
                               "Not analysed in GLASS" = "black",
                               "Different Results" = "red",
                               "Below Threshold" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Comparison with GLASS data") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

p
# data for table ####

# jarid2 genes

merge <- merge[!(is.na(merge$jarid2.gene)),]

nrow(merge)
table(merge$jarid2.gene)
table(merge$jarid2.tss)
table(merge$jarid2.tss == FALSE & merge$jarid2.gene == TRUE)

# all

r1.1 <- table(merge$padj < 0.05)[2]
r1.2 <- table(merge$padj < 0.05 & merge$log2FoldChange > 0)[2]
r1.3 <- table(merge$padj < 0.05 & merge$log2FoldChange < 0)[2]

r1.4 <- table(merge$threshold == TRUE)[2]
r1.5 <- table(merge$threshold == TRUE & merge$log2FoldChange > 0)[2]
r1.6 <- table(merge$threshold == TRUE & merge$log2FoldChange < 0)[2]

# non-jarid2

r2.1 <- table(merge$padj < 0.05 & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE)[2]
r2.2 <- table(merge$padj < 0.05 & merge$log2FoldChange > 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE)[2]
r2.3 <- table(merge$padj < 0.05 & merge$log2FoldChange < 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE)[2]

r2.4 <- table(merge$threshold == TRUE & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE)[2]
r2.5 <- table(merge$threshold == TRUE & merge$log2FoldChange > 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE)[2]
r2.6 <- table(merge$threshold == TRUE & merge$log2FoldChange < 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE)[2]


# jarid2 gene and tss

r3.1 <- table(merge$padj < 0.05 & merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE)[2]
r3.2 <- table(merge$padj < 0.05 & merge$log2FoldChange > 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE)[2]
r3.3 <- table(merge$padj < 0.05 & merge$log2FoldChange < 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE)[2]

r3.4 <- table(merge$threshold == TRUE & merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE)[2]
r3.5 <- table(merge$threshold == TRUE & merge$log2FoldChange > 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE)[2]
r3.6 <- table(merge$threshold == TRUE & merge$log2FoldChange < 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE)[2]

# jarid2 gene only

r4.1 <- table(merge$padj < 0.05 & merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE)[2]
r4.2 <- table(merge$padj < 0.05 & merge$log2FoldChange > 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE)[2]
r4.3 <- table(merge$padj < 0.05 & merge$log2FoldChange < 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE)[2]

r4.4 <- table(merge$threshold == TRUE & merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE)[2]
r4.5 <- table(merge$threshold == TRUE & merge$log2FoldChange > 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE)[2]
r4.6 <- table(merge$threshold == TRUE & merge$log2FoldChange < 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE)[2]


# jarid2 tss only

r5.1 <- table(merge$padj < 0.05 & merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE)[2]
r5.2 <- table(merge$padj < 0.05 & merge$log2FoldChange > 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE)[2]
r5.3 <- table(merge$padj < 0.05 & merge$log2FoldChange < 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE)[2]

r5.4 <- table(merge$threshold == TRUE & merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE)[2]
r5.5 <- table(merge$threshold == TRUE & merge$log2FoldChange > 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE)[2]
r5.6 <- table(merge$threshold == TRUE & merge$log2FoldChange < 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE)[2]

# make a table

test <- data.frame(Total1 = c(r1.1,r2.1,r3.1,r4.1,r5.1),
                   Upregulated1 = c(r1.2,r2.2,r3.2,r4.2,r5.2),
                   Downregulated2 = c(r1.3,r2.3,r3.3,r4.3,r5.3),
                   Total2 = c(r1.4,r2.4,r3.4,r4.4,r5.4),
                   Upregulated2 = c(r1.5,r2.5,r3.5,r4.5,r5.5),
                   Downregulated2 = c(r1.6,r2.6,r3.6,r4.6,r5.6))

test

write_clip(test)

cont <- data.frame("Not Significant" = c(1451, 379), "Significant" = c(2440, 1103))
cont

chisq.test(cont)


       
# comparing to gene-level results #####

# list of genes with significant DEIs

merge <- merge[!is.na(merge$Gene.EnsID.Simp),]

genes.sig.dei <- merge[merge$padj < 0.05,10]
genes.sig.dei <- unique(genes.sig.dei)

num.isoforms.sig <- merge[merge$padj < 0.05,]
num.isoforms.sig <- num.isoforms.sig %>% group_by(Gene.EnsID.Simp) %>% summarise(counted = n())
num.isoforms.sig$jarid2 <- num.isoforms.sig$Gene.EnsID.Simp %in% jarid2.gene


yes.jarid2 <- num.isoforms.sig[num.isoforms.sig$jarid2 == TRUE,]
hist(yes.jarid2$counted)
table(yes.jarid2$counted)
mean(yes.jarid2$counted)


no.jarid2 <- num.isoforms.sig[num.isoforms.sig$jarid2 == FALSE,]
hist(no.jarid2$counted)
table(no.jarid2$counted)
mean(no.jarid2$counted)


test <- as.data.frame(genes.sig.dei)
test$num.isoforms <- sum(merge$Gene.EnsID.Simp == test$genes.sig.dei, na.rm = TRUE)

merge$any.isoform.sig <- merge$Gene.EnsID.Simp %in% genes.sig.dei

length(unique(merge[merge$any.isoform.sig == TRUE & merge$gene.status == TRUE & merge$jarid2.gene == FALSE,10]))

length(unique(merge$Gene.EnsID.Simp))

head(test)

# MA-plot

# DESeq2::plotMA(res)
# shrunken LFC
# reslfc <- lfcShrink(dds, "tumourtype_R_vs_P", type = "apeglm")
# DESeq2::plotMA(resLFC)


# comparing with glass data #####

test <- merge[-log10(merge$padj.x) < 10 & abs(merge$log2FoldChange.x > 1),]
test <- test[!(is.na(test$padj.y)),]

plot(test$padj.x, test$padj.y)
plot(test$padj.x, test$padj.y, xlim = c(0,1), ylim = c(0,1))
hist(test$padj.y, breaks = c(0,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
cor(test$padj.x, test$padj.y)

head(test)

