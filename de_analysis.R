# import libraries

library(edgeR)

# import count data

dat <- read.delim("PvR_isoformCounts_all.txt", header = TRUE)

# check dimensions

dim(dat)                      

# rearrange columns so that isoform identifier, gene name and gene type in the first 3 columns, followed by the counts across 176 samples (for 88 patients)

dat <- dat[,c(1,10,11,2:9,12:ncol(dat))]
head(dat)

# remove rows that with counts summing to zero (i.e. no expression for transcript across all samples)

dat <- dat[rowSums(dat[,4:ncol(dat)]) != 0,]

# create DGEList using edgeR package

y <- DGEList(counts=dat[,4:179], genes=dat[,1:3])

# TMM normalisation

y <- calcNormFactors(y)
y$samples

# create design matrix 

# vector of all sample names in format patientid_tumourtype
samples <- rownames(y$samples)

# extract patient IDs and format as factors

patient <- sub("_.*","",samples)
patient <- factor(patient)

# extract tumour types (i.e. primary or recurrent) and format as factors

tumour.type <- sub(".*_","",samples)
tumour.type <- ifelse(tumour.type == "Primary","P", ifelse(tumour.type == "Recurrent", "R", tumour.type))
tumour.type <- as.factor(tumour.type)

# create design matrix

design <- model.matrix(~ patient + tumour.type)
rownames(design) <- colnames(y)
design

# estimate dispersions

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# do differential transcript expression

fit <- glmFit(y, design)
lrt <- glmLRT(fit)

# create dataframe of results and write to csv

res <- topTags(lrt)
write.table(res, "results.csv")


# total number of differentially expression transcripts

summary(de <- decideTestsDGE(lrt))

# plot results

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")


