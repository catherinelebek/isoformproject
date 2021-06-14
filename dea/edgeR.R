# import libraries

library(edgeR)

# import count data

dat <- read.delim("/nobackup/bs20chlb/inputdata/PvR_isoformCounts_all.txt", header = TRUE)

# rearrange columns

dat <- dat[,c(1,10,11,2:9,12:ncol(dat))]

# check columns in correct order

colnames(dat)

# create DGEList data class
# columns 4 to final column are counts per sample
# columns 1 to 3 are transcript ID, gene name and gene type

y <- DGEList(counts=dat[,4:ncol(dat)], genes=dat[,1:3])

# remove columns corresponding to samples that should be excluded from the analyses

# import list of patients to remove

patients.remove <- read.delim("/nobackup/bs20chlb/inputdata/patients_remove.txt", header = FALSE)

# convert list to vector

patients.remove <- as.vector(t(patients.remove))

# filter DGEList using the samples object

keep <- !sub("_.*","",rownames(y$samples)) %in% patients.remove # indexes for patient samples not in the patient remove file
y <- y[,keep] # subsetting y

# remove rows that add up to zero

keep <- rowSums(y$count) != 0 # indexes for transcripts with non-zero reads for at least one sample
y <- y[keep,] # subsetting y

# TMM normalisation

y <- calcNormFactors(y)

# create design matrix 

samples <- rownames(y$samples)
patient <- sub("_.*","",samples)
patient <- factor(patient)
tumour.type <- sub(".*_","",samples)
tumour.type <- ifelse(tumour.type == "Primary","P", ifelse(tumour.type == "Recurrent", "R", tumour.type))
tumour.type <- as.factor(tumour.type)
design <- model.matrix(~ patient + tumour.type)
rownames(design) <- colnames(y)

# estimate common dispersion

y <- estimateGLMCommonDisp(y, design)

# print out dispersion estimate to validate this stage in code was reached

cd <- y$common.dispersion

write(cd, "/nobackup/bs20chlb/outputdata/commondispersion.txt")

# estimate trended dispersion

y <- estimateGLMTrendedDisp(y, design)

# estimate tagwise dispersion

y <- estimateGLMTagwiseDisp(y, design)

# plot dispersions

plotBCV(y)

# do differential transcript expression

# fitting transcript-wise general linear models

fit <- glmFit(y, design)

# conduct likelihood ratio test for last coefficient in linear model (by default)

lrt <- glmLRT(fit)

# print results of top 10 most differentially expressed transcripts 

res <- topTags(lrt)
write.table(res, "/nobackup/bs20chlb/outputdata/edgeRresults.csv")

# total number of differentially expression transcripts

summary(de <- decideTestsDGE(lrt))

# plot results

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
