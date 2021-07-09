library(topGO)
library(biomaRt)
library(Rgraphviz)
library(clipr)

# connect to ensembl

bm <- useMart("ensembl")
bm <- useDataset("hsapiens_gene_ensembl", mart = bm)
EG2GO <- getBM(mart = bm, attribute = c("ensembl_gene_id", "ensembl_gene_id_version", "ensembl_transcript_id", "ensembl_transcript_id_version","external_gene_name", "go_id"))

# import results from isoform DEA

# isoforms <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv", 
#                         header = T, sep = ",")
isoforms <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv", 
                       header = T, sep = ",")
# isoforms <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results.csv", 
#                       header = T, sep = ",")

isoforms$Row.names.simp <- sub("\\..*","",isoforms$Row.names)

# get ensembl gene ids for transcripts in case need it for analysis in future
# mapping <- EG2GO[,c(1,3)]
# mapping <- unique(mapping)
# isoforms <- merge(isoforms, mapping, 
#                    by.x = "Row.names.simp", by.y = "ensembl_transcript_id", all.x = TRUE)


# analysis at the transcript level #####

# get the GO ID for each transcript

gid2GO <- by(EG2GO$go_id, EG2GO$ensembl_transcript_id, function(x) as.character(x))

# get adjusted p-value for each transcript

gL <- isoforms$padj

# label each adjusted p-value with ensembl transcript id

names(gL) <- isoforms$Row.names.simp
gL <- gL[!is.na(gL)]

# return transcripts with a p-value less than 0.05

topDiffGenes <- function(allScore) {return (allScore < 0.05)}
x <- topDiffGenes(gL)
table(x)

# create topGO data object
# allGenes are all the genes used in differential expression analysis
# geneSel are top differentiated genes
# gene2GO is a mapping between transcript identifiers and GO terms


GOdata <- new("topGOdata", ontology = "CC", allGenes=gL, geneSel = topDiffGenes,
              annot = annFUN.gene2GO, gene2GO = gid2GO, description = "All diff transcripts")


# Fisher test

resultFisher.classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultsFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

# KS test

resultsKS.classic <- runTest(GOdata, algorithm = "classic", statistic  = "ks")
resultsKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

# view results

allRes <- GenTable(GOdata, classicFisher=resultFisher.classic, classicKS=resultsKS.classic,
                   elimFisher = resultsFisher.elim, elimKS = resultsKS.elim, 
                   orderBy = "classicKS", ranksOf = "classicFisher", topNodes = 10)
allRes


write_clip(allRes)

# plot results

showSigOfNodes(GOdata, score(resultFisher.classic), firstSigNodes = 5, useInfo = "all")

# comparing test results ####

# function used in below plots to vary colours of points based on number of significantly differentially expressed genes

colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}

# KS classic vs. KS elim

pVal_KS_classic <- score(resultsKS.classic)
pVal_KS_elim <- score(resultsKS.elim)[names(pVal_KS_classic)]
gstat <- termStat(GOdata, names(pVal_KS_classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pVal_KS_classic, pVal_KS_elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

# Compare Fisher classic with Fisher elim
pVal_Fisher_classic <- score(resultFisher.classic)
pVal_Fisher_elim <- score(resultsFisher.elim)[names(pVal_Fisher_classic)]
gstat <- termStat(GOdata, names(pVal_Fisher_classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pVal_Fisher_classic, pVal_Fisher_elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)

# Compare Fisher classic with KS classic
pVal_Fisher_classic <- score(resultFisher.classic)
pVal_KS_classic <- score(resultsKS.classic)[names(pVal_Fisher_classic)]
gstat <- termStat(GOdata, names(pVal_Fisher_classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pVal_Fisher_classic, pVal_KS_classic, xlab = "p-value fisher classic", ylab = "p-value KS classic",
     pch = 19, cex = gSize, col = gCol)

# Compare Fisher elim with KS elim
pVal_Fisher_elim <- score(resultsFisher.elim)
pVal_KS_elim <- score(resultKS.elim)[names(pVal_Fisher_elim)]
gstat <- termStat(GOdata, names(pVal_Fisher_elim))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)
plot(pVal_Fisher_elim, pVal_KS_elim, xlab = "p-value classic", ylab = "p-value elim",
     pch = 19, cex = gSize, col = gCol)
