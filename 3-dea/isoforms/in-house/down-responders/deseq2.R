library(DESeq2)
# library(BiocParallel)

# counts <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/down-responders/PvR_isoformCounts_filtered.txt",header = T, sep = "\t")
counts <- read.delim("/nobackup/bs20chlb/inputdata/seconddata/down-responders/PvR_isoformCounts_filtered.txt",header = T, sep = "\t")

rownames(counts) <- counts[,1]
counts <- counts[,c(-1,-2,-3)]

samples <- data.frame(matrix(ncol = 2, nrow = ncol(counts)))
colnames(samples) <- c("patientid","tumourtype")
rownames(samples) <- colnames(counts)

for (i in 1:nrow(samples)){
  samples[i,1] <- gsub(".{2}$","",colnames(counts)[i])
  samples[i,2] <- substr(colnames(counts)[i], nchar(colnames(counts)[i]), nchar(colnames(counts)[i]))
}

samples$patientid <- as.factor(samples$patientid)
samples$tumourtype <- as.factor(samples$tumourtype)

all(rownames(samples) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = ~ patientid + tumourtype)

dds$tumourtype <- relevel(dds$tumourtype, ref = "P")

# dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(4))

dds <- DESeq(dds)

save(dds, file = "/nobackup/bs20chlb/outputdata/dea/seconddata/down-responders/deseq2.RData")
