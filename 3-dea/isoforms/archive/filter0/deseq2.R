library(DESeq2)

counts <- read.delim("/nobackup/bs20chlb/inputdata/filter0/PvR_isoformCounts_filtered.txt",header = T, sep = " ")
rownames(counts) <- counts[,1]
counts <- counts[,c(-1,-2,-3)]


samples <- data.frame(matrix(ncol = 2, nrow = ncol(counts)))
colnames(samples) <- c("patientid","tumourtype")
rownames(samples) <- colnames(counts)

for (i in 1:nrow(samples)){
  samples[i,1] <- sub("_.*","",colnames(counts)[i])
  samples[i,2] <- sub(".*_","",colnames(counts)[i])
}

samples$patientid <- as.factor(samples$patientid)
samples$tumourtype <- as.factor(samples$tumourtype)

all(rownames(samples) == colnames(counts))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = samples,
                              design = ~ patientid + tumourtype)

featureData <- data.frame(gene = counts[,2])
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)

dds$tumourtype <- relevel(dds$tumourtype, ref = "P")

dds <- DESeq(dds)

save(dds, file = "deseq.RData")

res <- results(dds)

resOrdered <- res[order(res$pvalue),]

summary(res)

sum(res$padj < 0.05, na.rm = TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)

