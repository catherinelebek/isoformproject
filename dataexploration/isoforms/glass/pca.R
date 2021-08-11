# import library for pca plots

library(factoextra)
library(ggsci)

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

# remove outliers

counts.delta <- counts.delta[,colnames(counts.delta) != "GLSS.MD.0023.R1" &
                               colnames(counts.delta) != "GLSS.MD.0022.R1" &
                               colnames(counts.delta) != "GLSS.SM.R083.R1"]

dim(counts.delta)
# transpose counts matrix

counts.delta.t <- t(counts.delta)
rownames(counts.delta.t)

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

pca.res <- prcomp(counts.delta.t, scale. = TRUE)


# plot results starting with PC1 vs. PC2

p <- fviz_pca_ind(pca.res, geom.ind = c("point", "text"),
                  axes = c(5,6),
                  pointsize = 3,  
                  title = " ",
                  habillage = responder.type$responder,
                  palette = c("#190F97", "#de3b04"),
                  repel = TRUE,
                  axes.linetype = NA,
                  mean.point = FALSE)


ggpubr::ggpar(p,
              font.x = c(15,"bold"),
              font.y = c(15,"bold"),
              font.legend = 15,
              font.tickslab = c(15),
              xlab = "PC3", ylab = "PC4",
              legend.title = "Responder-type", legend.position = "top",
              ggtheme = theme_bw()) +
              scale_shape_manual(values = c(16,16)) +           
              ggpubr::rremove("grid")

# Make a scree plot

fviz_screeplot(pca.res, addlabels = TRUE, ncp = 15,
               main = "Scree plot of the first 15 PCs",
               ggtheme = theme_minimal())


# extract isoforms

head(pca.res)
topisoforms <- pca.res$rotation
topisoforms <- as.data.frame(topisoforms)

# order the loadings for PC3

topisoforms <- topisoforms[order(topisoforms$PC3, decreasing = T),]

write.csv(topisoforms, "~/Documents/Semester3/Project/Results/resultsanalysis/PCA/topisoforms_GLASS.csv")

head(topisoforms)




