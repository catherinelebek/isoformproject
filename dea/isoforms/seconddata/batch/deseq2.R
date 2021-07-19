library(DESeq2)


counts <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/seconddata/isoform_BatchCorrected_PC_12072021.txt",
                     header = T, sep = "\t")
rownames(counts) <- counts$EnsID
counts <- counts[,c(-1)]
head(counts)

metadata <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/seconddata/MetaData_GT_250621.txt", 
                       sep = "\t", header = T)

metadata$ResponderType <- ifelse(metadata$NES < 0, "U","D")

samples <- data.frame(matrix(ncol = 2, nrow = ncol(counts)))
colnames(samples) <- c("patientid","tumourtype")
rownames(samples) <- colnames(counts)

for (i in 1:nrow(samples)){
  samples[i,1] <- gsub(".{2}$","",colnames(counts)[i])
  samples[i,2] <- substr(colnames(counts)[i], nchar(colnames(counts)[i]), nchar(colnames(counts)[i]))
}

samples <- merge(samples, metadata[,c("Patient.ID","LibraryType","Sample.Source","Location.Primary","PrimarySubtype",
                                     "RecurrentSubtype","ResponderType")], by.x = "patientid", by.y = "Patient.ID",
                                      all.x = TRUE)

samples <- samples[match(colnames(counts), paste0(samples$patientid,"_",samples$tumourtype)),]

rownames(samples) <- paste0(samples$patientid, "_", samples$tumourtype)
samples

all(rownames(samples) == colnames(counts))

# let's subset into primary and recurrent tumours
# expected 70 primary and 73 recurrent
# for this to work, the row order in samples must equal the column order in counts

idx.primary <- samples$tumourtype == "P"
samples.primary <- samples[idx.primary,]
counts.primary <- counts[,idx.primary]

idx.recurrent <- samples$tumourtype == "R"
samples.recurrent <- samples[idx.recurrent,]
counts.recurrent <- counts[,idx.recurrent]

# create subset using just Stranded_Total RNA

idx <- samples.primary$LibraryType == "Stranded_Total"

samples.primary.strandedtotal <- samples.primary[idx,]
counts.primary.strandedtotal <- counts.primary[,idx]

all(rownames(samples.primary.strandedtotal) == colnames(counts.primary.strandedtotal))

head(counts.primary.strandedtotal)

# create subset using just mRNA

idx <- grep("mRNA", samples.primary$LibraryType)

samples.primary.mRNAtotal <- samples.primary[idx,]
counts.primary.mRNAtotal <- counts.primary[,idx]

all(rownames(samples.primary.mRNAtotal) == colnames(counts.primary.mRNAtotal))

# create DESeq2 dataset

dds <- DESeqDataSetFromMatrix(countData = counts.primary.mRNAtotal,
                              colData = samples.primary.mRNAtotal,
                              design = ~ ResponderType)

# variance stabilisting transformed data in order to run PCA
dds.sf <- estimateSizeFactors(dds)
counts[1:6,1:3]
counts(dds)[1:6,1:3]
counts(dds.sf)[1:6,1:3]
counts(dds.sf, normalized = TRUE)[1:6,1:3]
boxplot(counts(dds.sf, normalized = TRUE))
vst <- varianceStabilizingTransformation(dds.sf)
boxplot(assay(vst))


# identifying variance in top most differentially expressed genes

plotPCA(vst, intgroup = c("ResponderType"))

# running DESeq2

dds$ResponderType <- relevel(dds$ResponderType, ref = "D")

dds <- DESeq(dds)

# import gene list for annotating results with gene symbols

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

# reorder

merge <- merge[order(merge$pvalue),]

# move gene list to be second column

merge <- merge[,c(1,8,2:7)]

# get a summary of the results at the 90% confidence level

summary(res)

# remove NAs

merge <- merge[!is.na(merge$padj),]

# create a threshold column

merge$threshold <- merge$padj < 0.05 & abs(merge$log2FoldChange) > 1

# volcano plot

ggplot(merge) +
  geom_point(aes(x = log2FoldChange, y=-log10(padj), colour = threshold)) +
  geom_text_repel(aes(x = log2FoldChange, y=-log10(padj),
                      label = ifelse(threshold == TRUE, GeneName, "")), size = 3) +
  ggtitle("Differential Isoform Expression - All Patients") +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  scale_colour_manual(values=c(1,2)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
  theme_bw()

head(merge)

write.csv(merge, "/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/deseq2results_mRNAtotal.csv")

