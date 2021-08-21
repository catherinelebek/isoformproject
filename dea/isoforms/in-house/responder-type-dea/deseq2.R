library(DESeq2)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(ggvenn)


counts <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/seconddata/isoform_BatchCorrected_PC_12072021.txt",
                     header = T, sep = "\t")
rownames(counts) <- counts$EnsID
counts <- counts[,c(-1)]
head(counts)

metadata <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/seconddata/MetaData_GT_250621.txt", 
                       sep = "\t", header = T)

metadata$ResponderType <- ifelse(metadata$NES < 0, "D","U")

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

table(samples$ResponderType, samples$tumourtype)

# create subset using just Stranded_Total RNA and without 3 outliers based on PC2

idx <- samples.primary$LibraryType == "Stranded_Total" &
       samples.primary$patientid != "Walton8" &
       samples.primary$patientid != "Walton10" &
       samples.primary$patientid != "Preston22"


samples.primary.strandedtotal.nooutliers <- samples.primary[idx,]
counts.primary.strandedtotal.nooutliers <- counts.primary[,idx]

all(rownames(samples.primary.strandedtotal) == colnames(counts.primary.strandedtotal))

head(counts.primary.strandedtotal)

# create subset using just mRNA

idx <- grep("mRNA", samples.primary$LibraryType)

samples.primary.mRNAtotal <- samples.primary[idx,]
counts.primary.mRNAtotal <- counts.primary[,idx]

all(rownames(samples.primary.mRNAtotal) == colnames(counts.primary.mRNAtotal))

# create subset using just DFKZ & Rabadan

idx <- samples.primary$Sample.Source == "DFKZ" | 
       samples.primary$Sample.Source == "Rabadan"

samples.primary.sources <- samples.primary[idx,]
counts.primary.sources <- counts.primary[,idx]

all(rownames(samples.primary.sources) == colnames(counts.primary.sources))


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

p <- plotPCA(vst, intgroup = c("Sample.Source"))
p + geom_text_repel(aes(label = name))

plotPCA(vst, intgroup = c("ResponderType"))
head(samples.primary)

plotPCA(vst, intgroup = c("Location.Primary"), pcs = c(3,4))

table(samples.primary.mRNAtotal$ResponderType)



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

# write to csv

write.csv(merge, paste0("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/deseq2results_mRNAtotal.csv"))

# create a threshold column

merge$threshold <- merge$padj < 0.05 & abs(merge$log2FoldChange) > 1

# volcano plot

p <- ggplot(merge) +
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

p

# results comparison

sources <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/deseq2results_sources.csv",
                    header = T, sep = ",")

all <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/deseq2results_all.csv",
                header = T, sep = ",")

mRNA <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/deseq2results_mRNAtotal.csv",
                 header = T, sep = ",")

strandedtotal <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/deseq2results_strandedtotal.csv",
                 header = T, sep = ",")



all_res <- all[,2:3]
colnames(all_res) <- c("EnsID","GeneName")

all_res <- merge(all_res, sources[,c("Row.names","log2FoldChange","padj")], by.x = "EnsID", by.y = "Row.names")
colnames(all_res)[(ncol(all_res)-1):ncol(all_res)] <- c("LFC.sources", "padj.sources")

all_res <- merge(all_res, all[,c("Row.names","log2FoldChange","padj")], by.x = "EnsID", by.y = "Row.names")
colnames(all_res)[(ncol(all_res)-1):ncol(all_res)] <- c("LFC.all", "padj.all")

all_res <- merge(all_res, mRNA[,c("Row.names","log2FoldChange","padj")], by.x = "EnsID", by.y = "Row.names")
colnames(all_res)[(ncol(all_res)-1):ncol(all_res)] <- c("LFC.mRNA", "padj.mRNA")

all_res <- merge(all_res, strandedtotal[,c("Row.names","log2FoldChange","padj")], by.x = "EnsID", by.y = "Row.names")
colnames(all_res)[(ncol(all_res)-1):ncol(all_res)] <- c("LFC.strandedtotal", "padj.strandedtotal")

head(all_res)

all_res[grep("COL5A2",all_res$GeneName),]

SOD2_transcripts
SOD2_transcripts <- all_res[grep("SOD2",all_res$GeneName),1:2]
SOD2_transcripts$EnsID <- sub("\\..","",SOD2_transcripts$EnsID)
write.table(SOD2_transcripts,
            "~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/SOD2.txt",
            row.names = F)



COL5A2_transcripts <- all_res[grep("COL5A2",all_res$GeneName),1:2]
COL5A2_transcripts$EnsID <- sub("\\..","",COL5A2_transcripts$EnsID)
write.table(COL5A2_transcripts,
            "~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/COL5A2.txt",
            row.names = F)

COL5A2_transcripts

NRXN3_transcripts <- all_res[grep("NRXN3",all_res$GeneName),1:2]
NRXN3_transcripts$EnsID <- sub("\\..","",NRXN3_transcripts$EnsID)
write.table(NRXN3_transcripts,
            "~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/NRXN3.txt",
            row.names = F)

FRYL_transcripts <- all_res[grep("FRYL",all_res$GeneName),1:2]
FRYL_transcripts$EnsID <- sub("\\..","",FRYL_transcripts$EnsID)
write.table(FRYL_transcripts,
            "~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/FRYL.txt",
            row.names = F)

FRYL_transcripts

# volcano plot

all_res$threshold <- ifelse(abs(all_res$LFC.sources) > 1 & all_res$padj.sources < 0.05,TRUE, FALSE)
all_res$direction <- ifelse(all_res$threshold == TRUE & all_res$LFC.sources > 1, "Up",
                            ifelse(all_res$threshold == TRUE & all_res$LFC.sources <  -1, "Down",
                                   "Below threshold"))

all_res$sources <- all_res$padj.sources < 0.05
all_res$all <- all_res$padj.all < 0.05
all_res$mRNA <- all_res$padj.mRNA < 0.05
all_res$strandedtotal <- all_res$padj.strandedtotal < 0.05
all_res$all.sig <- all_res$sources == TRUE & all_res$all == TRUE &
                    all_res$mRNA == TRUE & all_res$strandedtotal == TRUE
all_res$all.sources.stranded <- all_res$sources == TRUE & all_res$all == TRUE &
                                 all_res$strandedtotal == TRUE
all_res$all.sources <- all_res$sources == TRUE & all_res$all == TRUE



p <- ggplot(all_res, aes(x = LFC.sources, y = -log10(padj.sources))) +
  geom_point(aes(colour = direction, shape = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(all_res, all.sources.stranded == TRUE),
                   aes(label = GeneName, fill = direction),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  geom_text_repel(data = subset(all_res, -log10(padj.sources) > 9),
                   aes(label = GeneName),
                   min.segment.length = unit(0, "lines"),
                   segment.colour = "black") +
  scale_color_manual(values = c("Up" = "darkgreen",
                                "Down" = "#F11406", 
                                "Below threshold" = "grey")) +
  scale_fill_manual(values =  c("Up" = "darkgreen",
                                "Down" = "#F11406", 
                                "Below threshold" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Direction of dyregulation") +
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

ggplot(all_res, aes(A = mRNA, B = sources, C = strandedtotal, D = all)) +
  geom_venn(set_names = c("mRNA",
                          "DFKZ/Rabadan",
                          "Stranded_Total",
                          "All")) +
  theme_bw() +
  ggpubr::rremove("grid") +
  ggpubr::rremove("axis")



all_res$stranded.mRNA <- all_res$strandedtotal == TRUE & all_res$mRNA == TRUE
all_res$threshold <- ifelse(abs(all_res$LFC.strandedtotal) > 1 & all_res$padj.strandedtotal < 0.05,TRUE, FALSE)
all_res$direction <- ifelse(all_res$threshold == TRUE & all_res$LFC.strandedtotal > 1, "Up",
                            ifelse(all_res$threshold == TRUE & all_res$LFC.strandedtotal <  -1, "Down",
                                   "Below threshold"))


p <- ggplot(all_res, aes(x = LFC.strandedtotal, y = -log10(padj.strandedtotal))) +
  geom_point(aes(colour = direction, shape = threshold)) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p-value") +
  geom_vline(xintercept = 1, linetype = 2) +
  geom_vline(xintercept = -1, linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_label_repel(data = subset(all_res, stranded.mRNA == TRUE),
                   aes(label = GeneName, fill = direction),
                   fontface = "bold", color = "white",
                   min.segment.length = unit(0, "lines"),
                   box.padding = unit(0.1, "lines"),
                   segment.colour = "black") +
  geom_text_repel(data = subset(all_res, -log10(padj.strandedtotal) > 9),
                  aes(label = GeneName),
                  min.segment.length = unit(0, "lines"),
                  segment.colour = "black") +
  scale_color_manual(values = c("Up" = "darkgreen",
                                "Down" = "#F11406", 
                                "Below threshold" = "grey")) +
  scale_fill_manual(values =  c("Up" = "darkgreen",
                                "Down" = "#F11406", 
                                "Below threshold" = "grey")) +
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19)) +
  guides(shape = "none", fill = "none") +
  labs(color = "Direction of dyregulation") +
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
