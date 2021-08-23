library(DESeq2)
# library(BiocParallel)

counts <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/PvR_isoformCounts_filtered.txt",header = T, sep = "\t")
# counts <- read.delim("/nobackup/bs20chlb/inputdata/seconddata/PvR_isoformCounts_filtered.txt",header = T, sep = "\t")
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

p <- plotPCA(vst, intgroup = "Sample.Source")
p
pca.res <- p$data

# comparing PC1 coordinates

pca.pc1 <- data.frame(matrix(ncol = 2, nrow = nrow(pca.res)/2))
patients <- unique(sub(".{2}$","",rownames(pca.res)))
rownames(pca.pc1) <- patients
colnames(pca.pc1) <- c("P","R")

for (i in 1:length(patients)){
  pca.pc1[i,1] <- pca.res[pca.res$patientid == patients[i] &
                          pca.res$tumourtype == "P",1]
  pca.pc1[i,2] <- pca.res[pca.res$patientid == patients[i] &
                            pca.res$tumourtype == "R",1]
}

plot(pca.pc1$P, pca.pc1$R, pch = 16, col = "darkblue", xlab = "Primary PC1", ylab = "Recurrent PC1")
 cor(pca.pc1$P, pca.pc1$R)
t.test(pca.pc1$P, pca.pc1$R, paired = TRUE, alternative = "two.sided", mu = 0, 
       var.equal = TRUE, conf.level= 0.95)

# comparing PC2 coordinates

pca.pc2 <- data.frame(matrix(ncol = 2, nrow = nrow(pca.res)/2))
patients <- unique(sub(".{2}$","",rownames(pca.res)))
rownames(pca.pc2) <- patients
colnames(pca.pc2) <- c("P","R")

for (i in 1:length(patients)){
  pca.pc2[i,1] <- pca.res[pca.res$patientid == patients[i] &
                            pca.res$tumourtype == "P",2]
  pca.pc2[i,2] <- pca.res[pca.res$patientid == patients[i] &
                            pca.res$tumourtype == "R",2]
}

plot(pca.pc2$P, pca.pc2$R, pch = 16, col = "darkblue", xlab = "Primary PC2", ylab = "Recurrent PC2")
cor(pca.pc2$P, pca.pc2$R)
t.test(pca.pc1$P, pca.pc1$R, paired = TRUE, alternative = "two.sided", mu = 0, 
       var.equal = TRUE, conf.level= 0.95)

# running DESeq2

dds$tumourtype <- relevel(dds$tumourtype, ref = "P")

dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(4))



save(dds, file = "/nobackup/bs20chlb/outputdata/dea/seconddata/deseq2.RData")
