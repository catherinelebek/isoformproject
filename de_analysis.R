# import libraries
# checking git setup correctly 

library(edgeR)

# import count data

dat <- read.delim("PvR_isoformCounts_all.txt", header = TRUE)
dim(dat)                      

# rearrange columns

dat <- dat[,c(1,10,11,2:9,12:ncol(dat))]
head(dat)

# remove rows that add up to zero

dat <- dat[rowSums(dat[,4:ncol(dat)]) != 0,]

# create DGEList data class

y <- DGEList(counts=dat[,4:179], genes=dat[,1:3])

# TMM normalisation

y <- calcNormFactors(y)
y$samples

# create design matrix 

samples <- rownames(y$samples)
patient <- sub("_.*","",samples)
patient <- factor(patient)
tumour.type <- sub(".*_","",samples)
tumour.type <- ifelse(tumour.type == "Primary","P", ifelse(tumour.type == "Recurrent", "R", tumour.type))
tumour.type <- as.factor(tumour.type)
data.frame(Samples = colnames(y), patient, tumour.type)
design <- model.matrix(~ patient + tumour.type)
rownames(design) <- colnames(y)
design

# estimate dispersion

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# do differential transcript expression

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
res <- topTags(lrt)
write.table(res, "results.csv")


# total number of differentially expression transcripts

summary(de <- decideTestsDGE(lrt))

# plot results

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")


