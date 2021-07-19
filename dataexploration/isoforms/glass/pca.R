# import library for pca plots

library(factoextra)
library(ggsci)

# import glass data

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
jarid

# create dataframe of sample and their corresponding responder type

samples <- as.data.frame(sub(".{3}$","",rownames(counts.delta.t)))
colnames(samples) <- "X"
samples

responder.type <- merge(samples, jarid[,c("Patient","responder.type","direction")], 
                        by.x = "X", by.y = "Patient")

responder.type <- responder.type[match(sub(".{3}$","",rownames(counts.delta.t)),responder.type$X),]
table(responder.type$X == sub(".{3}$","",rownames(counts.delta.t)))
table(responder.type$responder)
responder.type$responder <- as.factor(responder.type$responder)

# run pca

counts.delta.t <- t(counts.delta)
rownames(counts.delta.t)

pca.res <- prcomp(counts.delta.t)


# plot results starting with PC1 vs. PC2

p <- fviz_pca_ind(pca.res, geom = "point",
                  axes = c(1,2),
                  pointsize = 3,  
                  habillage = responder.type$responder.type,
                  repel = TRUE,
                  palette = "npg")

ggpubr::ggpar(p,
              title = "Principal Component Analysis",
              subtitle = "LFC Isoforms",
              caption = "Source: GLASS data",
              legend.title = "Responder Type", legend.position = "top")

# Make a scree plot

fviz_screeplot(pca.res, addlabels = TRUE, ncp = 15,
               main = "Scree plot of the first 15 PCs",
               ggtheme = theme_minimal())


# Hierarchical clustering

d <- dist(counts.delta.t)
h <- hclust(d, labels = c(t(responder.type$responder.type)))
h$labels <- responder.type$responder.type
plot(h)

responder.type
h$labels
responder.type$responder.type
