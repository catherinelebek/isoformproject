library(DESeq2)
library(ggplot2)
library(ggrepel)
library(biomaRt)
library(dplyr)

# load DESeq2 object as dds

load("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2.RData")

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

merge <- merge[order(res$pvalue),]

# move gene list to be second column

merge <- merge[,c(1,8,2:7)]

# get a summary of the results at the 90% confidence level

summary(res)

# remove NAs

merge <- merge[!is.na(merge$padj),]

# save results to csv

# write.csv(merge, "/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv")

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


# create volcano plot

merge$threshold <- merge$padj < 0.05 & abs(merge$log2FoldChange) > 1
merge$label <- ifelse(merge$threshold == 1 & merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE,
                      "Sig - JARID2 Gene & TSS",
                      ifelse(merge$threshold == 1 & merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE,
                      "Sig - JARID2 Gene Only",
                      ifelse(merge$threshold == 1 & merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE,
                      "Sig - JARID2 Isoform Only",
                      ifelse(merge$threshold == 1 & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE, 
                      "Sig",
                      ifelse(merge$threshold == 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == TRUE,
                      "Not Sig - JARID2 Gene & TSS",
                       ifelse(merge$threshold == 0 & merge$jarid2.gene == TRUE & merge$jarid2.tss == FALSE,
                      "Not Sig - JARID2 Gene Only",
                       ifelse(merge$threshold == 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == TRUE, 
                      "Not Sig - JARID2 Isoform Only",
                       ifelse(merge$threshold == 0 & merge$jarid2.gene == FALSE & merge$jarid2.tss == FALSE,
                      "Not Sig","null"))))))))
                                    

merge$top <- ifelse(-log10(merge$padj) > 5 | abs(merge$log2FoldChange) > 10,1,0)

merge$label <- as.factor(merge$label)

ggplot(merge) +
  geom_point(aes(x = log2FoldChange, y=-log10(padj), colour = label)) +
  geom_text_repel(aes(x = log2FoldChange, y=-log10(padj),
                label = ifelse(top == 1, GeneName, "")), size = 3) +
  ggtitle("Differential Isoform Expression - All Patients") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_colour_manual(values=c(8,7,6,5,4,3,2,1)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme_bw()

# plot specific genes

gene <- "EIF4E3" # set gene of interest

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


