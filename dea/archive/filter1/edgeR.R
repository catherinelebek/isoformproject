# import libraries

library(edgeR)

# import count data

dat <- read.table("/nobackup/bs20chlb/inputdata/nonstringent/PvR_isoformCounts_filtered.txt", header = TRUE)

# create DGEList data class

y <- DGEList(counts=dat[,4:ncol(dat)], genes=dat[,1:3])

# TMM normalisation

y <- calcNormFactors(y)

# create design matrix 

# create a vector with all the sample names - these are in the format PATIENTID_TUMOURTYPE

samples <- rownames(y$samples)

# extract patient ids from sample names and format as factors

patient <- sub("_.*","",samples)
patient <- as.factor(patient)

# extract tumour type from sample names, convert all to "P" or "R" format, and format as factors

tumour.type <- sub(".*_","",samples)
tumour.type <- as.factor(tumour.type)

# finally create the design matrix and rename the rows with the sample names

design <- model.matrix(~ patient + tumour.type)
rownames(design) <- colnames(y)

# estimate common dispersion

y <- estimateGLMCommonDisp(y, design)

# print out dispersion estimate to validate this stage in code was reached

cd <- y$common.dispersion

cd

# estimate trended dispersion

y <- estimateGLMTrendedDisp(y, design)

# estimate tagwise dispersion

y <- estimateGLMTagwiseDisp(y, design)

save(y, file = "tagwisedisp.RData")

# plot dispersions

plotBCV(y)

# do differential transcript expression

# fitting transcript-wise general linear models

fit <- glmFit(y, design)

# conduct likelihood ratio test for last coefficient in linear model (by default)

lrt <- glmLRT(fit)

# print results of top 10 most differentially expressed transcripts 

res <- topTags(lrt, n = Inf)
write.table(res, "/nobackup/bs20chlb/outputdata/dea/nonstringent/edgeRresults.csv")

# total number of differentially expression transcripts

summary(de <- decideTestsDGE(lrt))

# plot results

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
