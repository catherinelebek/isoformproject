library(DESeq2)
library(ggplot2)
library(ggrepel)
library(bioma)

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

# save results to csv

# write.csv(merge, "/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv")

# MA-plot

# DESeq2::plotMA(res)

# shrunken LFC
# reslfc <- lfcShrink(dds, "tumourtype_R_vs_P", type = "apeglm")
# DESeq2::plotMA(resLFC)

# import summarydf

bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart = bm)
EG2GO <- getBM(mart = bm, attribute = c("ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"))

jarid2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",")

head(jarid2)
EG2GO$JARID2 <- EG2GO$ensembl_gene_id %in% jarid2$V1

head(EG2GO)

EG2GO[EG2GO$external_gene_name == "SCN8A",]

# ggplot

merge$EnsIDSimp <- sub("\\..*","",merge$Row.names)

merge <- merge(merge, EG2GO[,c("ensembl_transcript_id", "JARID2")], 
                      by.x = "EnsIDSimp", by.y = "ensembl_transcript_id",
                      all.x = TRUE)



merge$threshold <- merge$padj < 0.05 & abs(merge$log2FoldChange) > 1
merge$threshold <- ifelse(merge$threshold == 1 & merge$JARID2 == TRUE, "Sig - JARID2", 
                          ifelse(merge$threshold == 1, "Sig", "Not Sig"))

merge <- merge[!is.na(merge$padj),]
merge <- merge[!is.na(merge$threshold),]

merge$top <- ifelse(-log10(merge$padj) > 5 | abs(merge$log2FoldChange) > 10,1,0)

merge$threshold <- as.factor(merge$threshold)

ggplot(merge) +
  geom_point(aes(x = log2FoldChange, y=-log10(padj), colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, y=-log10(padj),
                label = ifelse(top == 1, GeneName, "")), size = 3) +
  ggtitle("Differential Isoform Expression - All Patients") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_colour_manual(values=c("black","red", "green")) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme_bw()

