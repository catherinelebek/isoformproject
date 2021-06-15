# import libraries

library(edgeR)

# import count data

dat <- read.delim("/nobackup/bs20chlb/inputdata/PvR_isoformCounts_all.txt", header = TRUE)

# rearrange columns

dat <- dat[,c(1,10,11,2:9,12:ncol(dat))]

# create DGEList data class

y <- DGEList(counts=dat[,4:ncol(dat)], genes=dat[,1:3])

# remove columns corresponding to samples that should be excluded from the analyses

# import list of patients to remove

patients.remove <- read.delim("/nobackup/bs20chlb/inputdata/patients_remove.txt", header = F)

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

# create a vector with all the sample names - these are in the format PATIENTID_TUMOURTYPE

samples <- rownames(y$samples)

# extract patient ids from sample names and format as factors

patient <- sub("_.*","",samples)
patient <- factor(patient)

# extract tumour type from sample names, convert all to "P" or "R" format, and format as factors

tumour.type <- sub(".*_","",samples)
tumour.type <- ifelse(tumour.type == "Primary","P", ifelse(tumour.type == "Recurrent", "R", tumour.type))
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

restop <- topTable(fit, coef = "tumour.typeR")

# write results to a file

write.table(restop, "/nobackup/bs20chlb/outputdata/limmaresults.csv")

plotMDS(y, labels = patient, col = as.numeric(tumour.type))

