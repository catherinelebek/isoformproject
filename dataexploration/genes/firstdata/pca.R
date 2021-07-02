library(stringr)

counts.full <- read.delim("~/Documents/Semester3/Project/InputData/genes/PvR_geneCounts_all_LS_23062021.txt.txt",header = T, sep = "\t")
head(counts.full)

new.patients <- read.csv("~/Documents/Semester3/Project/InputData/isoforms/firstdata/samplefilters/newpatients.csv", header = F)
new.patients <- c(t(new.patients))

counts <- counts.full

rownames(counts) <- counts[,1]
counts <- counts[,c(-1,-2,-3)]


colnames(counts) <- sub("Primary","P",colnames(counts))
colnames(counts) <- sub("Recurrent","R",colnames(counts))

colnames(counts)

remove <- !gsub(".{2}$","",colnames(counts)) %in% new.patients
counts <- counts[,remove]

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



counts.primary <- apply(counts.primary, c(1,2), log2)
counts.recurrent <- apply(counts.recurrent, c(1,2), log2)

counts.delta <- counts.recurrent - counts.primary

head(counts.delta)


dim(counts.delta)

test <- rowSums(abs(counts.delta)) == Inf
table(is.na(test))

# need patients to be the rows

counts.delta <- na.omit(counts.delta)
idx <- rowSums(abs(counts.delta)) != Inf
counts.delta <- counts.delta[idx,]


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

plot(pca.res$rotation[,1],pca.res$rotation[,2],xlab="loading 1",ylab="loading 2", 
     main="Loadings - LFC gene expression, 42 patients", col = pca.res.rot$type, pch = 16)
legend(0.2,0.4, c("Down","Up"), col = c(1,2), pch = 16)

