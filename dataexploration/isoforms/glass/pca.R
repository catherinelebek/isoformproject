library(stringr)

counts.full <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/isoforms/glass/PvR_isoformCounts_filtered.txt",
                          header = T, sep = "\t")
head(counts.full)

counts <- counts.full

rownames(counts) <- counts[,1]
counts <- counts[,c(-1,-2)]

func <- function(x){x + 0.0001}

counts <- apply(counts, 1:2, func)

head(counts)

cols <- gsub(".{3}$","",colnames(counts))
cols <- unique(cols)
cols

primary <- gsub(".TP$","",colnames(counts)) %in% cols
recurrent <- gsub(".R1$","",colnames(counts)) %in% cols

counts.primary <- counts[,primary]
counts.recurrent <- counts[,recurrent]

table(gsub(".TP$","",colnames(counts.primary)) == gsub(".R1$","",colnames(counts.recurrent)))

head(counts.primary)
head(counts.recurrent)

counts.primary <- apply(counts.primary, c(1,2), log2)
counts.recurrent <- apply(counts.recurrent, c(1,2), log2)
counts.delta <- counts.recurrent - counts.primary

head(counts.delta)


dim(counts.delta)

table(rowSums(abs(counts.delta)) == Inf)
table(is.na(counts.delta))

# counts.delta <- na.omit(counts.delta)
# idx <- rowSums(abs(counts.delta)) != Inf
# counts.delta <- counts.delta[idx,]

# import jarid2 response type data

jarid <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/glassdata/JARID2_results.txt",
                     header = TRUE, sep = "\t")

# get responder type status

jarid$NES <- as.numeric(jarid$NES)
jarid$direction <- ifelse(is.na(jarid$NES), jarid$ES, jarid$NES)
jarid$responder.type <- ifelse(jarid$direction > 0, "Up", "Down")
table(jarid$responder.type)
head(jarid)

# run pca

pca.res <- prcomp(counts.delta, scale. = T)
pca.res.rot <- pca.res$rotation
pca.res.rot <- as.data.frame(pca.res.rot)
rownames(pca.res.rot) <- gsub(".{3}$","",rownames(pca.res.rot))

# merge data


pca.res.rot <- merge(pca.res.rot, jarid[,c("Patient","responder.type")], by.x = "row.names", by.y = "Patient", all.x = TRUE)
pca.res.rot$responder.type <- as.factor(pca.res.rot$responder.type)

plot(pca.res.rot$PC1,pca.res.rot$PC2,xlab="loading 1",ylab="loading 2", 
     main="Loadings - LFC isoform expression", col = pca.res.rot$responder.type, pch = 16)
legend(-0.4,0.4,c("Down","Up"), col = c(1,2), pch = 16)
