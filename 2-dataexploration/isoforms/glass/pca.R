# This script runs a PCA on the log2FC in normalised isoform counts between 66 paired samples from the glass dataset

# import library for pca plots

library(factoextra)
library(ggsci)
library(patchwork)

# import glass data

counts.full <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/isoforms/glass/glassfilter/PvR_isoformCounts_filtered.txt",
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

# transpose counts matrix

counts.delta.t <- t(counts.delta)
rownames(counts.delta.t)

# create dataframe of sample and their corresponding responder type

samples <- as.data.frame(sub(".{3}$","",rownames(counts.delta.t)))
colnames(samples) <- "X"
samples

responder.type <- merge(samples, jarid[,c("Patient","direction")], 
                        by.x = "X", by.y = "Patient")

responder.type <- responder.type[match(sub(".{3}$","",rownames(counts.delta.t)),responder.type$X),]
table(responder.type$X == sub(".{3}$","",rownames(counts.delta.t)))
responder.type



# run pca

pca.res <- prcomp(counts.delta.t, scale. = TRUE)


# plot results with PC3 vs. PC4

p <- fviz_pca_ind(pca.res, geom = "point",
                  axes = c(3,4),
                  title = " ",
                  pointsize = 3,  
                  col.ind = responder.type$direction,
                  repel = TRUE,
                  gradient.cols = c("#3933FF", "#E7B800", "#FC4E07"),
                  axes.linetype = NA)

p1 <- ggpubr::ggpar(p,
              font.x = c(15,"bold"),
              font.y = c(15,"bold"),
              font.legend = 15,
              font.tickslab = c(15),
              xlab = "PC3", ylab = "PC4",
              legend.title = "NES", legend.position = "top",
              ggtheme = theme_bw()) +
  ggpubr::rremove("grid")

p1 


# Make a scree plot


p <- fviz_screeplot(pca.res, 
                    addlabels = TRUE, 
                    ncp = 10,
                    title = " ",
                    barfill = "#190F9799",
                    barcolor = "#190F97",
                    xlab = "PC")

p2 <- ggpubr::ggpar(p,
              font.x = c(15, "bold"),
              font.y = c(15, "bold"),
              font.tickslab = c(12),
              ggtheme = theme_bw()) +
  ggpubr::rremove("grid")

p2

p1 + p2

# extract isoforms

head(pca.res)
topisoforms <- pca.res$rotation
topisoforms <- as.data.frame(topisoforms)

# order the loadings for PC3

topisoforms <- topisoforms[order(topisoforms$PC3, decreasing = T),]

write.csv(topisoforms, "~/Documents/Semester3/Project/Results/resultsanalysis/PCA/topisoforms_GLASS.csv")

head(topisoforms)




