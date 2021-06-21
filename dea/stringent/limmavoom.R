# import libraries

library(edgeR)

# import count data

dat <- read.table("/nobackup/bs20chlb/inputdata/stringent/PvR_isoformCounts_filtered.txt", header = TRUE)

# create DGEList data class

y <- DGEList(counts=dat[,4:ncol(dat)], genes=dat[,1:3])

# TMM normalisation

y <- calcNormFactors(y)

# create design matrix 

# create a vector with all the sample names - these are in the format PATIENTID_TUMOURTYPE

samples <- rownames(y$samples)

# extract patient ids from sample names and format as factors

patient <- sub("_.*","",samples)
patient <- factor(patient)

# extract tumour type from sample names, convert all to "P" or "R" format, and format as factors

tumour.type <- sub(".*_","",samples)
tumour.type <- as.factor(tumour.type)

# finally create the design matrix and rename the rows with the sample names

design <- model.matrix(~ patient + tumour.type)
rownames(design) <- colnames(y)

# run with limma

# voom converts read counts to log2(CPM) with associated weights
# creates an EList object

v <- voom(y, design, plot=TRUE)

# fit a linear model

fit <- lmFit(v, design)

# test for differentially expression transcripts

fit <- eBayes(fit)

# print out table of top most differentially expression transcripts

restop <- topTable(fit, coef = "tumour.typeR", n = Inf)

# write results to a file

write.table(restop, "/nobackup/bs20chlb/outputdata/dea/stringent/limmaresults.csv")

plotMDS(y, labels = patient, col = as.numeric(tumour.type))

write.table(restop, "/Users/catherinehogg/Documents/Semester3/Project/Results/dea/limmaresults.csv")


