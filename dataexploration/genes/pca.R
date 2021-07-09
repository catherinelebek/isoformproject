library(stringr)

counts.full <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/genes/PvR_geneCounts_filtered.txt",header = T, sep = "\t")

counts <- counts.full

rownames(counts) <- counts[,1]
counts <- counts[,c(-1,-2,-3)]
colnames(counts)

# run following lines if want to just do PCA using JARID2 genes ####

jarid2 <- read.delim("~/Documents/Semester3/Project/InputData/genes/jarid2.csv",header = F, sep = ",")
jarid2 <- c(t(jarid2))

counts <- counts[sub("\\..*","",rownames(counts)) %in% jarid2,]

######

func <- function(x){x + 0.01}

counts <- apply(counts, 1:2, func)

head(counts)

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

counts.primary <- apply(counts.primary, c(1,2), log2)
counts.recurrent <- apply(counts.recurrent, c(1,2), log2)

counts.delta <- counts.recurrent - counts.primary

head(counts.delta)


dim(counts.delta)

# check 

table(rowSums(abs(counts.delta)) == Inf)
table(is.na(counts.delta))

# need patients to be the rows
# counts.delta <- na.omit(counts.delta)
# idx <- rowSums(abs(counts.delta)) != Inf
# counts.delta <- counts.delta[idx,]


pca.res <- prcomp(counts.delta, scale. = T)
pca.res.rot <- pca.res$rotation
pca.res.rot <- as.data.frame(pca.res.rot)
rownames(pca.res.rot) <- gsub(".{2}$","",rownames(pca.res.rot))

metadata <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/MetaData_GT_250621.txt",
                                 header = TRUE, "\t")

head(metadata)

pca.res.rot <- merge(pca.res.rot, metadata[,c("Patient.ID","NES")], by.x = "row.names", by.y = "Patient.ID", all.x = TRUE)
pca.res.rot$type <- ifelse(pca.res.rot$NES > 0, "Up", "Down")
pca.res.rot$type <- as.factor(pca.res.rot$type)

plot(pca.res.rot$PC1,pca.res.rot$PC2,xlab="loading 1",ylab="loading 2", 
     main="Loadings - LFC gene expression", col = pca.res.rot$type, pch = 16)
legend(-0.2,-0.25,c("Down","Up"), col = c(1,2), pch = 16)


# find the most important genes to PC1

topgenes <- as.data.frame(pca.res$x)
topgenes <- as.data.frame(topgenes$PC1)
rownames(topgenes) <- rownames(pca.res$x)
topgenes <- apply(topgenes,1,abs)
topgenes <- topgenes[order(topgenes, decreasing = T)]
topgenes <- names(topgenes)

write.csv(topgenes, "/Users/catherinehogg/Documents/Semester3/Project/Results/resultsanalysis/topgenes.csv", row.names = F)


