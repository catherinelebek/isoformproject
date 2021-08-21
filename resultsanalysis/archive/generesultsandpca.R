# comparing signficicant results for all genes with those expected from JARID2 dysregulation #####


genes.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/deseq2results.csv",
                        header = T, sep = ",")

genes.all <- genes.all[,c(-1)]
colnames(genes.all)[1] <- "EnsID"

jarid2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                        header = F, sep = ",")
jarid2 <- c(t(jarid2))

# 3717 jarid2 genes tested for DEA

jarid2tested <- jarid2[jarid2 %in% sub("\\..*","",genes.all$EnsID)]

genes.all.sig <- genes.all[genes.all$padj <= 0.05,]
genes.all.not.sig <- genes.all[genes.all$padj > 0.05,]

genes.all.sig$jarid2 <- sub("\\..*","",genes.all.sig$EnsID) %in% jarid2
genes.all.not.sig$jarid2 <- sub("\\..*","",genes.all.not.sig$EnsID) %in% jarid2

# 383 jarid2 genes identified as being significantly dysregulated

table(genes.all.sig$jarid2)

# 3324 jarid2 genes not identified as being significantly dysregulated

table(genes.all.not.sig$jarid2)

# which of the significant genes are important in PC1?

topgenes <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/resultsanalysis/topgenes.csv")
topgenes <- c(t(topgenes))

head(topgenes)
genes.all.sig$topgenes <- match(genes.all.sig$EnsID, topgenes)
genes.all.not.sig$topgenes <- match(genes.all.not.sig$EnsID, topgenes)
table(genes.all.sig$topgenes)
head(genes.all.sig)

boxplot(genes.all.sig$topgenes)

# comparing signficicant results for down-responder genes with those expected from JARID2 dysregulation #####

genes.down.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/down-responders/deseq2results.csv",
                        header = T, sep = ",")

head(genes.down.all)

genes.down.all <- genes.down.all[,c(-1)]
colnames(genes.down.all)[1] <- "EnsID"

jarid2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",")
jarid2 <- c(t(jarid2))

# 3721 jarid2 genes tested for DEA

jarid2tested <- jarid2[jarid2 %in% sub("\\..*","",genes.up.all$EnsID)]

genes.down.all.sig <- genes.down.all[genes.down.all$padj <= 0.05,]
genes.down.all.not.sig <- genes.down.all[genes.down.all$padj > 0.05,]

genes.down.all.sig$jarid2 <- sub("\\..*","",genes.down.all.sig$EnsID) %in% jarid2
genes.down.all.not.sig$jarid2 <- sub("\\..*","",genes.down.all.not.sig$EnsID) %in% jarid2

# 782 jarid2 genes identified as being significantly dysregulated

table(genes.down.all.sig$jarid2)

# 3913 jarid2 genes not identified as being significantly dysregulated

table(genes.down.all.not.sig$jarid2)

# which of the significant genes are important in PC1?

topgenes <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/resultsanalysis/topgenes.csv")
topgenes <- c(t(topgenes))

head(topgenes)
genes.down.all.sig$topgenes <- match(genes.down.all.sig$EnsID, topgenes)
table(genes.down.all.sig$topgenes)
head(genes.down.all.sig, 100)



boxplot(genes.down.all.sig$topgenes)

# comparing signficicant results for up-responder genes with those expected from JARID2 dysregulation #####

genes.up.all <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/dea/genes/up-responders/deseq2results.csv",
                           header = T, sep = ",")

head(genes.up.all)

genes.up.all <- genes.up.all[,c(-1)]
colnames(genes.up.all)[1] <- "EnsID"

jarid2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",")
jarid2 <- c(t(jarid2))

# 3717 jarid2 genes tested for DEA

jarid2tested <- jarid2[jarid2 %in% sub("\\..*","",genes.up.all$EnsID)]

genes.up.all.sig <- genes.up.all[genes.up.all$padj <= 0.05,]
genes.up.all.not.sig <- genes.up.all[genes.up.all$padj > 0.05,]

genes.up.all.sig$jarid2 <- sub("\\..*","",genes.up.all.sig$EnsID) %in% jarid2
genes.up.all.not.sig$jarid2 <- sub("\\..*","",genes.up.all.not.sig$EnsID) %in% jarid2

# 1484 jarid2 genes identified as being significantly dysregulated

table(genes.up.all.sig$jarid2)

# 3913 jarid2 genes not identified as being significantly dysregulated

table(genes.up.all.not.sig$jarid2)

# which of the significant genes are important in PC1?

topgenes <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/resultsanalysis/topgenes.csv")
topgenes <- c(t(topgenes))

head(topgenes)
head(genes.up.all)
genes.up.all$topgenes <- match(genes.up.all$EnsID, topgenes)
table(genes.up.all$topgenes)
head(genes.up.all, 100)

table(is.na(genes.up.all.sig$topgenes))

plot(genes.up.all$padj, genes.up.all$topgenes)
cor(genes.up.all$padj, genes.up.all$topgenes)

genes.up.all.sig$topgenes <- match(genes.up.all.sig$EnsID, topgenes)
genes.up.all.not.sig$topgenes <- match(genes.up.all.not.sig$EnsID, topgenes)
genes.down.all.sig$topgenes <- match(genes.down.all.sig$EnsID, topgenes)
genes.down.all.not.sig$topgenes <- match(genes.down.all.not.sig$EnsID, topgenes)

par(mfrow = c(3,2))
boxplot(genes.all.sig$topgenes, main = "PC1-rankings for signficiant DE genes", cex.main = 0.7)
boxplot(genes.all.not.sig$topgenes, main = "PC1-rankings for not signficiant DE genes", cex.main = 0.7)
boxplot(genes.up.all.sig$topgenes, main = "PC1-rankings for significant DE genes in up-responders", cex.main = 0.7)
boxplot(genes.up.all.not.sig$topgenes, main = "PC1-rankings for not signficiant DE genes in up-responders", cex.main = 0.7)
boxplot(genes.down.all.sig$topgenes, main = "PC1-rankings for significant DE genes in down-responders", cex.main = 0.7)
boxplot(genes.down.all.not.sig$topgenes, main = "PC1-rankings for not signficiant DE genes in down-responders", cex.main = 0.7)



head(genes.all.sig)
head(genes.all.not.sig)


head(genes.down.all.sig)    