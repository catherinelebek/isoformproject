# import libraries

library(edgeR)

# setting file paths

intputdata <- file.path("/nobackup/bs20chlb","inputdata")
outputdata <- file.path("/nobackup/bs20chlb","outputdata")

# import count data

setwd(intputdata)
dat <- read.delim("PvR_isoformCounts_all.txt", header = TRUE)

# rearrange columns

dat <- dat[,c(1,10,11,2:9,12:ncol(dat))]

# create DGEList data class

y <- DGEList(counts=dat[,4:ncol(dat)], genes=dat[,1:3])

# remove columns corresponding to samples that should be excluded from the analyses

# import list of patients to remove

patients.remove <- read.delim("patients_remove.txt", header = FALSE)

# convert list to vector

patients.remove <- as.vector(t(patients.remove))

# filter DGEList using the samples object

keep <- !sub("_.*","",rownames(y$samples)) %in% patients.remove
y <- y[,keep]

# remove rows that add up to zero

keep <- rowSums(y$count) != 0
y <- y[keep,]

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

# estimate dispersion

y <- estimateGLMCommonDisp(y, design)

cd <- y$common.dispersion

setwd(outputdata)
write(cd, "commondispersion.txt")

y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# do differential transcript expression

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
res <- topTags(lrt)
write.table(res, "edgeRresults.csv")

# total number of differentially expression transcripts

summary(de <- decideTestsDGE(lrt))

# plot results

detags <- rownames(y)[as.logical(de)]
plotSmear(lrt, de.tags=detags)
abline(h=c(-1, 1), col="blue")
