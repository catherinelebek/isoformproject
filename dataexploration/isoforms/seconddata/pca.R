library(stringr)

counts.full <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/PvR_isoformCounts_filtered.txt",header = T, sep = "\t")
head(counts)

counts <- counts.full

rownames(counts) <- counts[,1]
counts <- counts[,c(-1,-2,-3)]

idx_primary <- colnames(counts)

colnames(counts)

cols <- gsub(".{2}$","",colnames(counts))
cols <- unique(cols)
cols

primary <- gsub("_P$","",colnames(counts)) %in% cols
recurrent <- gsub("_R$","",colnames(counts)) %in% cols

counts.primary <- counts[,primary]
counts.recurrent <- counts[,recurrent]

table(gsub("_P$","",colnames(counts.primary)) == gsub("_R$","",colnames(counts.recurrent)))

head(counts.primary)
head(counts.recurrent)

counts.primary <- is.


counts.primary.log <- apply(counts.primary, c(1,2), log2)
counts.recurrent.log <- apply(counts.recurrent, c(1,2), log2)
counts.delta.log <- counts.recurrent.log - counts.primary.log

counts.primary[1:5,1]
counts.recurrent[1:5,1]
counts.primary.log[1:5,1]
counts.recurrent.log[1:5,1]
counts.delta.log[1:5,1]


counts.delta <- counts.recurrent - counts.primary

counts.delta.log <- counts.delta.log[is.infinite(counts.delta.log)] <- -99999

head(counts.delta.log)

# need patients to be the rows


counts.delta <- na.omit(counts.delta)
idx <- rowSums(abs(counts.delta)) != Inf
counts.delta <- counts.delta[idx,]


pca.res <- prcomp(counts.delta)
pca.res.rot <- pca.res$rotation
pca.res.rot <- as.data.frame(pca.res.rot)
rownames(pca.res.rot) <- gsub(".{2}$","",rownames(pca.res.rot))

metadata <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/MetaData_GT_250621.txt",
                                 header = TRUE, "\t")

head(metadata)

pca.res.rot <- merge(pca.res.rot, metadata[,c("Patient.ID","NES")], by.x = "row.names", by.y = "Patient.ID", all.x = TRUE)
pca.res.rot$type <- ifelse(pca.res.rot$NES > 0, "Up", "Down")
pca.res.rot$type <- as.factor(pca.res.rot$type)

plot(pca.res$rotation[,1],pca.res$rotation[,2],xlab="loading 1",ylab="loading 2", 
     main="Loadings", col = pca.res.rot$type, pch = 16)

