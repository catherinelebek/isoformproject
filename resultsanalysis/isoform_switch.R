# for up-responders ####

degoppdir <- read.csv("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/degoppdir.csv", header = T)
head(degoppdir)


gene.counts <- gene.counts[gene.counts$Gene.EnsID %in% degoppdir$Gene.EnsID,]
gene.counts <- merge(gene.counts, degoppdir[,c("Gene.EnsID","GeneName")], by.x = "Gene.EnsID", by.y = "Gene.EnsID", all.x = TRUE)

gene.counts <- gene.counts[-1,]

gene.counts <- gene.counts[,c(21,1,3:20)]
gene.counts

isoform.switch <- gene.counts[gene.counts$threshold.up > 1,]
isoform.switch

isoform.switch[,c(1,8,9,13,14,17,18)]

# for down responders ####

degoppdir <- read.csv("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/degoppdir.csv", header = T)
head(degoppdir)

gene.counts <- read.csv("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/genecounts.csv", header = T)


gene.counts <- gene.counts[gene.counts$Gene.EnsID %in% degoppdir$Gene.EnsID,]
gene.counts <- merge(gene.counts, degoppdir[,c("Gene.EnsID","GeneName")], by.x = "Gene.EnsID", by.y = "Gene.EnsID", all.x = TRUE)

gene.counts <- gene.counts[,c(13,1,3:12)]
gene.counts

isoform.switch <- gene.counts[gene.counts$threshold.up > 1,]
isoform.switch
