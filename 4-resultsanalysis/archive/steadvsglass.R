stead.up <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv", header = T, sep = ",")
stead.down <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results.csv", header = T, sep = ",")
glass.up <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/glass/up-responders/deseq2results.csv", header = T, sep = ",")
glass.down <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/glass/down-responders/deseq2results.csv", header = T, sep = ",")

stead.up <- stead.up[,c(-1)]
stead.down <- stead.down[,c(-1)]
stead.up$EnsID <- sub("\\..*","",stead.up$Row.names)
stead.down$EnsID <- sub("\\..*","",stead.down$Row.names)

# check isoforms unique

test <- stead.down[,9]
length(unique(test))

# get list of all isoforms tested in up responders (whether in stead and/or glass data) ####

up.isoforms <- c(t(stead.up[,9]), t(glass.up[,1]))
up.isoforms <- unique(up.isoforms)
up.compare <- as.data.frame(matrix(ncol = 1, nrow = length(up.isoforms)))
up.compare[,1] <- up.isoforms

head(up.compare)

up.compare <- merge(up.compare, stead.up[,c("EnsID","log2FoldChange","pvalue","padj")], 
                    by.x = "V1", by.y = "EnsID", all.x = TRUE)

up.compare <- merge(up.compare, glass.up[,c("X","log2FoldChange","pvalue","padj")],
                    by.x = "V1", by.y = "X", all.x = TRUE)


colnames(up.compare) <- c("EnsID", "Stead.LFC", "Stead.pvalue", "Stead.padj", "GLASS.LFC", "GLASS.pvalue", "GLASS.padj")

head(up.compare)

up.compare$stead.up <- ifelse(up.compare$Stead.LFC > 0, T, F)
up.compare$glass.up <- ifelse(up.compare$GLASS.LFC > 0, T, F)
up.compare$stead.sig <- ifelse(up.compare$Stead.padj < 0.05, T, F)
up.compare$glass.sig <- ifelse(up.compare$GLASS.padj < 0.05, T, F)


test <- up.compare[up.compare$stead.sig == T & up.compare$glass.sig == T & !is.na(up.compare$GLASS.padj),]
test <-  up.compare[up.compare$stead.sig == T & up.compare$glass.sig == F & !is.na(up.compare$GLASS.LFC),]
test <-  up.compare[up.compare$stead.sig == T & is.na(up.compare$GLASS.LFC),]
test <-  up.compare[up.compare$glass.sig == T & !is.na(up.compare$glass.sig) & up.compare$stead.sig == F & !is.na(up.compare$Stead.LFC),]

table(test$stead.up, test$glass.up)

# get list of all isoforms tested in down responders (whether in stead and/or glass data) ####

down.isoforms <- c(t(stead.down[,9]), t(glass.down[,1]))
down.isoforms <- unique(down.isoforms)
down.compare <- as.data.frame(matrix(ncol = 1, nrow = length(down.isoforms)))
down.compare[,1] <- down.isoforms

head(down.compare)

down.compare <- merge(down.compare, stead.down[,c("EnsID","log2FoldChange","pvalue","padj")], 
                    by.x = "V1", by.y = "EnsID", all.x = TRUE)

down.compare <- merge(down.compare, glass.down[,c("X","log2FoldChange","pvalue","padj")],
                    by.x = "V1", by.y = "X", all.x = TRUE)


colnames(down.compare) <- c("EnsID", "Stead.LFC", "Stead.pvalue", "Stead.padj", "GLASS.LFC", "GLASS.pvalue", "GLASS.padj")

head(down.compare)


down.compare$stead.up <- ifelse(down.compare$Stead.LFC > 0, T, F)
down.compare$glass.up <- ifelse(down.compare$GLASS.LFC > 0, T, F)
down.compare$stead.sig <- ifelse(down.compare$Stead.padj < 0.05, T, F)
down.compare$glass.sig <- ifelse(down.compare$GLASS.padj < 0.05, T, F)


test <- down.compare[down.compare$stead.sig == T & down.compare$glass.sig == T & !is.na(down.compare$GLASS.padj),]
test <-  down.compare[down.compare$stead.sig == T & down.compare$glass.sig == F & !is.na(down.compare$GLASS.LFC),]
test <-  down.compare[down.compare$stead.sig == T & is.na(down.compare$GLASS.LFC),]
test <-  down.compare[down.compare$glass.sig == T & !is.na(down.compare$glass.sig) & down.compare$stead.sig == F & !is.na(down.compare$Stead.LFC),]

table(test$stead.up, test$glass.up)

# what are the isoforms that overlap? ####

# up-responders

test1 <- up.compare[up.compare$stead.sig == T & up.compare$glass.sig == T & !is.na(up.compare$GLASS.padj),]
stead.up[stead.up$EnsID %in% test$EnsID,]

# down-responders

test2 <- down.compare[down.compare$stead.sig == T & down.compare$glass.sig == T & !is.na(down.compare$GLASS.padj),]

stead.down[stead.down$EnsID %in% test$EnsID,]

# in up and down responders

test1[test1$EnsID %in% test2$EnsID,]

stead.up[stead.up$EnsID == "ENST00000264790",]
stead.down[stead.down$EnsID == "ENST00000264790",]

stead.up[stead.up$EnsID == "ENST00000366574",]
stead.down[stead.down$EnsID == "ENST00000366574",]
      