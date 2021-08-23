edger <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/filter3/edgeRresultsfilter3.csv", header = T, sep = ",")
deseq <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/filter3/deseq2results.csv", header = T, sep = ",")

head(edger)
head(deseq)

merge1 <- merge(deseqsig, edger, by.x = "Row.names", by.y = "EnsID", all.x = TRUE)
merge1

head(deseqsig)
edgersig <- edger[edger$FDR < 0.05,]
deseqsig <- deseq[!is.na(deseq$padj) & deseq$padj < 0.05, ]

test <- unique(c(t(sub("-.*","",deseqsig$GeneName))))
test <- paste(test, collapse = ", ")
test

merge <- merge(deseqsig, edgersig, by.x = "Row.names", by.y = "EnsID")

merge <- merge[,c(1,3,5,13,8,16,9,17)]
merge


genes <- merge$GeneName.x
genes
genes <- paste(genes, collapse=", ")
genes

deseq["EIF4E3" == sub("-.*","",deseq$GeneName),]

