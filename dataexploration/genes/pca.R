# import library for pca plots

library(factoextra)

# import transcript filtered isoform counts data for all patients

counts.full <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/genes/PvR_geneCounts_filtered.txt",
                          header = T, sep = "\t")
head(counts.full)

counts <- counts.full

# make rownames the transcript EnsIDs

rownames(counts) <- counts[,1]

# remove the EnsIDs, GeneType and GeneType columns

counts <- counts[,c(-1,-2,-3)]

# function to add 0.01 to each value

func <- function(x){x + 0.01}

# apply function

counts <- apply(counts, 1:2, func)

head(counts)

# remove the tumour type from each patient id and collate in new vector, cols

cols <- gsub(".{2}$","",colnames(counts))
cols <- unique(cols)
cols

# split data into primary and recurrent tumour samples

primary <- gsub("_P$","",colnames(counts)) %in% cols
recurrent <- gsub("_R$","",colnames(counts)) %in% cols

counts.primary <- counts[,primary]
counts.recurrent <- counts[,recurrent]

# check the ordering of columns in each sample is the same

table(gsub("_P$","",colnames(counts.primary)) == gsub("_R$","",colnames(counts.recurrent)))

# take the log2 of each count value and then calculated the LFC

counts.primary <- apply(counts.primary, c(1,2), log2)
counts.recurrent <- apply(counts.recurrent, c(1,2), log2)
counts.delta <- counts.recurrent - counts.primary

# Make sure there are no NA and inifinite values

table(rowSums(abs(counts.delta)) == Inf)
table(is.na(counts.delta))

# transpose dataset into wide format for PCA analysis 

counts.delta.t <- t(counts.delta)

# remove tumour type from rownames

rownames(counts.delta.t) <- sub(".{2}$","",rownames(counts.delta.t))
dim(counts.delta.t)

# run PCA

pca.res <- prcomp(counts.delta.t, scale. = F)

# import metadata

metadata <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/MetaData_GT_250621.txt",
                     header = TRUE, "\t")

# create dataframe of sample and their corresponder NES scores

samples <- as.data.frame(rownames(counts.delta.t))
colnames(samples) <- "X"

responder.type <- merge(samples, metadata[,c("Patient.ID","NES")], 
                        by.x = "X", by.y = "Patient.ID")

responder.type <- responder.type[match(rownames(counts.delta.t),responder.type$X),]
table(responder.type$X == rownames(counts.delta.t))
table(responder.type$X == rownames(pca.res$x))
responder.type$responder <- ifelse(responder.type$NES < 0, "Down", "Up")
table(responder.type$responder)
responder.type$responder <- as.factor(responder.type$responder)

# plot results starting with PC1 vs. PC2

p <- fviz_pca_ind(pca.res, geom = "point",
                  axes = c(1,2),
                  pointsize = 3,  
                  col.ind = responder.type$NES,
                  repel = TRUE,
                  gradient.cols = c("#3933FF", "#E7B800", "#FC4E07"))

ggpubr::ggpar(p,
              title = "Principal Component Analysis",
              subtitle = "LFC Genes",
              caption = "Source: Stead Data",
              xlab = "PC1", ylab = "PC2",
              legend.title = "NES", legend.position = "top")


# then plot PC3 vs PC4

p <- fviz_pca_ind(pca.res, geom = "point",
                  axes = c(3,4),
                  pointsize = 3,  
                  col.ind = responder.type$NES,
                  repel = TRUE,
                  gradient.cols = c("#3933FF", "#E7B800", "#FC4E07"))

ggpubr::ggpar(p,
              title = "Principal Component Analysis",
              subtitle = "LFC Genes",
              caption = "Source: Stead Data",
              xlab = "PC3", ylab = "PC4",
              legend.title = "NES", legend.position = "top")


# Make a scree plot

fviz_screeplot(pca.res, addlabels = TRUE, ncp = 15,
               main = "Scree plot of the first 15 PCs",
               ggtheme = theme_minimal())



# extract most important genes to PC2

head(pca.res)
topgenes <- pca.res$rotation
topgenes <- as.data.frame(topgenes)
topgenes$PC2 <- abs(topgenes$PC2)
topgenes <- topgenes[order(topgenes[,2], decreasing = TRUE),]
write.csv(topgenes, "/Users/catherinehogg/Documents/Semester3/Project/Results/resultsanalysis/topgenes.csv",
          row.names = T)


