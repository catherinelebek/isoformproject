library(topGO)
library(biomaRt)
library(Rgraphviz)
library(clipr)

responder.type <- "up-responders"
direction <- "up"

# connect to ensembl

bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart = bm)
EG2GO <- getBM(mart = bm, attribute = c("ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version","external_gene_name", "go_id"))

# import results from isoform DEA

dea.res <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/",responder.type,"/deseq2results.csv"), 
                         header = T, sep = ",")

head(dea.res)
dea.res$Row.names <- sub("\\..*","",dea.res$Row.names)
dea.res$interest <- ifelse(dea.res$padj < 0.05 & dea.res$log2FoldChange > 1, "Sig", "Not Sig")

gL <- c(t(dea.res$interest))
gL <- as.factor(gL)
names(gL) <- dea.res$Row.names

head(gL)

# return transcripts with a p-value less than 0.05

topDiffGenes <- function(allScore) {return (allScore == "Sig")}
x <- topDiffGenes(gL)
table(x)

# get ensembl gene ids for transcripts in case need it for analysis in future
# mapping <- EG2GO[,c(1,3)]
# mapping <- unique(mapping)
# isoforms <- merge(isoforms, mapping, 
#                    by.x = "Row.names.simp", by.y = "ensembl_transcript_id", all.x = TRUE)


# analysis at the transcript level #####

# get the GO ID for each transcript

gid2GO <- by(EG2GO$go_id, EG2GO$ensembl_transcript_id, function(x) as.character(x))



# create topGO data object
# allGenes are all the genes used in differential expression analysis
# geneSel are top differentiated genes
# gene2GO is a mapping between transcript identifiers and GO terms


GOdata <- new("topGOdata", ontology = "BP", allGenes=gL, geneSel = topDiffGenes,
              annot = annFUN.gene2GO, gene2GO = gid2GO, description = "All diff transcripts",
              nodeSize = 5)

# Fisher test

resultFisher.classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultsFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# KS test

resultsKS.classic <- runTest(GOdata, algorithm = "classic", statistic  = "ks")
resultsKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
resultsKS.weight <- runTest(GOdata, algorithm = "weight01", statistic = "ks")

# view results

allRes <- GenTable(GOdata, classicFisher=resultFisher.classic, weightKS=resultsKS.weight, 
                   orderBy = "classicFisher", ranksOf = "weightKS", topNodes = 500)

allRes <- GenTable(GOdata, classicKS=resultsKS.classic, elimKS = resultsKS.elim,
                   orderBy = "classicKS", topNodes = 1000)


head(allRes, 20)
allRes[grep("RNA", allRes$Term),]

write_clip(allRes)

# plot results

showSigOfNodes(GOdata, score(resultFisher.classic), firstSigNodes = 20, useInfo = "all")