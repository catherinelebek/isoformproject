# DESeq2 first data vs DESeq2 second data ###################

res_isoforms_1 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/filter3/firstdata/deseq2results.csv", header = T, sep = ",")
res_isoforms_2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/filter3/seconddata/deseq2results.csv", header = T, sep = ",")

res_isoforms_1 <- res_isoforms_1[,-1]
res_isoforms_2 <- res_isoforms_2[,-1]

rownames(res_isoforms_1) <- res_isoforms_1[,1]
rownames(res_isoforms_2) <- res_isoforms_2[,1]

res_isoforms_1 <- res_isoforms_1[,-1]
res_isoforms_2 <- res_isoforms_2[,-1]

head(res_isoforms_1)
head(res_isoforms_2)

merge <- merge(res_isoforms_1, res_isoforms_2, by = "row.names", all.x = TRUE)

merge <- merge[order(merge$padj.x),]

merge$Both <- ifelse(merge$padj.x < 0.05 & merge$padj.y < 0.05, 1,0)
both <- merge[merge$Both == 1,]
both <- both[!is.na(both$Both),]


merge$First <- ifelse(merge$padj.x < 0.05 & (merge$padj.y >= 0.05 | is.na(merge$padj.y)),1,0)
first <- merge[merge$First == 1,]
first <- both[!is.na(both$Both),]


# genes data

res_genes <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/deseq2results.csv", header = T, sep = ",")


head(res_genes)


res_genes[res_genes$GeneName == "SCN8A",]


# glass data

res_glass <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/glass/deseq2results.csv", header = T, sep = ",")
gene_list_full <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt", header = T, sep = "\t")

gene_list <- gene_list_full

gene_list <- gene_list[,c(1,10)]
gene_list$EnsID <- gsub("\\..*.","",gene_list$EnsID)

head(res_glass)

res_glass <- merge(res_glass, gene_list, by.x = "X", by.y = "EnsID", all.x = TRUE)

res_glass <- res_glass[order(res_glass$padj),]
head(res_glass,100)
