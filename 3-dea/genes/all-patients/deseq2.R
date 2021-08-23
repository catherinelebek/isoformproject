# This script runs gene-level DEA using DESeq2 on the paired samples from 66 patients (in-house dataset)


library(DESeq2)
library(BiocParallel)
library(clipr)

counts <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/genes/PvR_geneCounts_filtered.txt",header = T, sep = "\t")
# counts <- read.delim("/nobackup/bs20chlb/inputdata/seconddata/filter3/PvR_isoformCounts_filtered.txt",header = T, sep = "\t")
rownames(counts) <- counts[,1]
counts <- counts[,c(-1,-2,-3)]

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

samples <- merge(samples, metadata[,c("Patient.ID","Sample.Source","Location.Primary","Overall.Survival..months.",
                                      "ResponderType","PrimarySubtype","RecurrentSubtype",
                                       "LibraryType")],
                 by.x = "patientid", by.y = "Patient.ID", all.x = TRUE)

rownames(samples)

samples <- samples[match(colnames(counts), paste0(samples$patientid,"_",samples$tumourtype)),]
rownames(samples) <- paste0(samples$patientid,"_",samples$tumourtype)
head(samples)

samples$patientid <- as.factor(samples$patientid)
samples$tumourtype <- as.factor(samples$tumourtype)
samples$Sample.Source <- as.factor(samples$Sample.Source)
samples$Location.Primary <- as.factor(samples$Location.Primary)
samples$ResponderType <- as.factor(samples$ResponderType)
samples$PrimarySubtype <- as.factor(samples$PrimarySubtype)
samples$RecurrentSubtype <- as.factor(samples$RecurrentSubtype)
samples$LibraryType <- as.factor(samples$LibraryType)


all(rownames(samples) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = ~ patientid + tumourtype)

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

plotPCA(vst, intgroup = c("LibraryType"))
write_clip(p)
    
dds$tumourtype <- relevel(dds$tumourtype, ref = "P")

dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(4))

save(dds, file = "~/Documents/Semester3/Project/Results/dea/genes/deseq2.RData")

table(samples$LibraryType)

# split dataset

dds1 <- dds[,samples$LibraryType == "Stranded_Total"]
dds2 <- dds[,samples$LibraryType != "Stranded_Total"]

dds1.sf <- estimateSizeFactors(dds1)
dds1.vst <- varianceStabilizingTransformation(dds1.sf)

dds2.sf <- estimateSizeFactors(dds2)
dds2.vst <- varianceStabilizingTransformation(dds2.sf)


# identifying variance in top most differentially expressed genes

plotPCA(dds1.vst, intgroup = c("ResponderType"))

plotPCA(dds2.vst, intgroup = c("ResponderType"))


# hierarchical clustering

d <- dist(t(assay(dds1.vst)))
h <- hclust(d)
plot(h)

# k-means clustering

k <- kmeans(t(assay(dds1.vst)), centers = 2)
k$cluster

plot(k$cluster)


