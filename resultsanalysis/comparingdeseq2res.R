resfilter0 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/filter0/deseq2results.csv", header = T, sep = ",")
resfilter3 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/filter3/deseq2results.csv", header = T, sep = ",")

resfilter0sig <- resfilter0[!is.na(resfilter0$padj) & resfilter0$padj < 0.05,]
resfilter3sig <- resfilter3[!is.na(resfilter3$padj) & resfilter3$padj < 0.05,]


head(resfilter0sig)

write.csv(resfilter0sig, "/Users/catherinehogg/Documents/Semester3/Project/Results/resultsanalysis/resfilter0sig.csv")
write.csv(resfilter3sig, "/Users/catherinehogg/Documents/Semester3/Project/Results/resultsanalysis/resfilter3sig.csv")


resfilter3[resfilter3$Row.names == "ENST00000547430.2",]

missing <- resfilter0sig[!resfilter0sig$Row.names %in% resfilter3sig$Row.names,2]

missingfilter3 <- resfilter3[resfilter3$Row.names %in% missing,]
missingfilter3
