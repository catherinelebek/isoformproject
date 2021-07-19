library(DESeq2)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(dplyr)
library(ggsci)

# load DESeq2 object as dds

load("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2.RData")

# load full list of transcripts in order to pull through gene names

genelist <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

# extract just transcripts and gene names

genelist <- genelist[,c(1,10)]

# create results table from dds object

res <- results(dds, alpha = 0.05)

# order results by raw p-value

resOrdered <- res[order(res$pvalue),]

# put in a dataframe so can merge with the list of transcripts and gene names

resOrdered <- as.data.frame(resOrdered)

# run merge to pull through gene names into the DESeq2 results table

merge <- merge(resOrdered, genelist, by.x = "row.names", by.y = "EnsID", all.x = TRUE)

# again order by p-value

merge <- merge[order(merge$pvalue),]

# move gene list to be second column

merge <- merge[,c(1,8,2:7)]

# get a summary of the results at the 90% confidence level

summary(res)

# remove NAs

merge <- merge[!is.na(merge$padj),]

# save results to csv

write.csv(merge, "/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv")

# MA-plot

# DESeq2::plotMA(res)
# shrunken LFC
# reslfc <- lfcShrink(dds, "tumourtype_R_vs_P", type = "apeglm")
# DESeq2::plotMA(resLFC)

# annotate transcripts with whether from a JARID2 gene

# import jarid2 gene list
jarid2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",")

# use ensembl database to attribute an ensembl transcript name to jarid2 gene Ens ID
# this is required so that the results from DESeq2 can be matched with JARID2 gene status
# using Ens transcript ID as an identifier

bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart = bm)
EG2GO <- getBM(mart = bm, attribute = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"))

jarid2 <- merge(jarid2, EG2GO[,c("ensembl_gene_id","external_gene_name")], by.x = "V1",
                by.y = "ensembl_gene_id", all.x = TRUE)

jarid2 <- jarid2 %>% distinct(external_gene_name)

merge$GeneIDSimp <- sub("-.*","",merge$GeneName)
merge$jarid2.gene <- merge$GeneIDSimp %in% jarid2$external_gene_name


# now need to label each transcript from DESeq2 results with information on whether it has a JARID2 start site
# import jarid 2 tss transcripts

jarid2.tss <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/jarid.tss.transcripts.csv",
                       header = T)

jarid2.tss <- c(t(jarid2.tss))

merge$jarid2.tss <- merge$Row.names %in% jarid2.tss

head(merge)

# create volcano plot

merge$threshold <- merge$padj < 0.05 & abs(merge$log2FoldChange) > 1
merge$jarid.status <- ifelse(merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE, "Gene & TSS",
                             ifelse(merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE, "Gene",
                                    ifelse(merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE, "TSS","Neither")))

merge$top <- ifelse(-log10(merge$padj) > 5 | abs(merge$log2FoldChange) > 10,1,0)

merge$jarid.status <- as.factor(merge$jarid.status)
merge$threshold <- as.factor(merge$threshold)

p <- ggplot(merge) +
  geom_point(aes(x = log2FoldChange, y=-log10(padj), colour = threshold, shape = jarid.status)) +
  geom_text_repel(aes(x = log2FoldChange, y=-log10(padj),
                      label = ifelse(top == 1, GeneName, "")), size = 3) +
  ggtitle("Differential Isoform Expression - Up-Responders") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_color_manual(values = c("#E7B800", "#FC4E07")) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme_bw()

p

# plot specific genes #####

gene <- "SPOCK3" # set gene of interest

merge.gene <- merge[grep(gene, merge$GeneName),]
merge.gene$top <- 1
merge.gene


ggplot(merge.gene) +
  geom_point(aes(x = log2FoldChange, y=-log10(padj), colour = label)) +
  geom_text_repel(aes(x = log2FoldChange, y=-log10(padj),
                      label = ifelse(top == 1, GeneName, "")), size = 3) +
  ggtitle("Differential Isoform Expression - All Patients") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_colour_manual(values=c(1,2,3,4,5,6,7)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme_bw()


# plot jarid only genes ######


merge.jaridonly <- merge[merge$jarid2.gene == TRUE | merge$jarid2.tss == TRUE,]


ggplot(merge.jaridonly) +
  geom_point(aes(x = log2FoldChange, y=-log10(padj), colour = label)) +
  geom_text_repel(aes(x = log2FoldChange, y=-log10(padj),
                      label = ifelse(top == 1, GeneName, "")), size = 3) +
  ggtitle("Differential Isoform Expression - Up-Responders") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_colour_manual(values=c(1,2,3,4,5,6,7)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# mean padj and LFC for JARID2 and non-JARID2 #####

padj1 <- merge[merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE,]
table(padj1$threshold)
padj1 <- padj1[padj1$threshold == TRUE,]
padj.mean1 <- mean(padj1$padj)
lfc.mean1 <- mean(padj1$log2FoldChange)
padj.mean1
lfc.mean1

padj2 <- merge[merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE,]
table(padj2$threshold)
padj2 <- padj2[padj2$threshold == TRUE,]
padj.mean2 <- mean(padj2$padj)
lfc.mean2 <- mean(padj2$log2FoldChange)
padj.mean2
lfc.mean2


t.test(padj1$padj, padj2$padj, alternative = c("two.sided"), mu = 0, 
       var.equal = FALSE, conf.level = 0.95, paired = FALSE)

t.test(padj1$log2FoldChange, padj2$log2FoldChange, alternative = c("two.sided"), mu = 0, 
       var.equal = FALSE, conf.level = 0.95, paired = FALSE)


# chi-squared test

test <- merge[merge$jarid2.tss == FALSE,]
cont.table <- table(test$threshold, test$jarid2.gene == TRUE)
cont.table 
chisq.test(cont.table)

# running tests to for association with direction of dysregulation

merge$direction <- ifelse(merge$log2FoldChange < 0, "Down", "Up")

nrow(merge[merge$threshold == TRUE & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE & merge$direction == "Down",])


# compare to glass results ######

head(merge)

glass.up <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/glass/up-responders/deseq2results.csv",
                       sep = ",", header = T)

merge$EnsIDSimp <- sub("\\..*","",merge$Row.names)
head(merge)

merge <- merge(merge, glass.up[,c("X","padj","log2FoldChange")], by.x = "EnsIDSimp", by.y = "X", all.x = TRUE)
merge$glass.up <- ifelse(merge$padj.y < 0.05 & merge$log2FoldChange.x*merge$log2FoldChange.y > 0, "Same", 
                         ifelse(merge$padj.y < 0.05, "Sig not same", FALSE))
merge$glass.up <- ifelse(is.na(merge$glass.up), "Not analysed in GLASS", merge$glass.up)

ggplot(merge) +
  geom_point(aes(x = log2FoldChange.x, y=-log10(padj.x), colour = glass.up)) +
  geom_text_repel(aes(x = log2FoldChange.x, y=-log10(padj.x),
                      label = ifelse(glass.up == TRUE, GeneName, "")), size = 3) +
  ggtitle("Differential Isoform Expression - Up-Responders") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_colour_manual(values=c(1,2,3,4)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme_bw()

merge[merge$glass.up == "Sig not same",]

