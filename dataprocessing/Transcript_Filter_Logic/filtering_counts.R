# import libraries

library(edgeR)

datfull.counts <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/PvR_isoformCounts_all.txt", header = TRUE)

dat <- datfull.counts

# import metadata

metadata <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/samplefilters/Metadata.csv", header = TRUE)

# import list of patients to remove based on metadata values

patients.remove <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/samplesfilters/patients_remove.txt", header = FALSE)

# convert to vector

patients.remove <- as.vector(t(patients.remove))

# import list of patients to remove based on reads < 30m

below30 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/samplefilters/below30.txt", header = FALSE)

# covert datamframe to vector

below30 <- as.vector(t(below30))

# rearrange columns

c0 <- grep("EnsID", colnames(dat))
c1 <- grep("GeneName", colnames(dat))
c2 <- grep("GeneType", colnames(dat))
c3 <- c1 - 1
c4 <- c2 + 1
dat <- dat[,c(c0,c1,c2,2:c3,c4:ncol(dat))]
colnames(dat)

dim(dat)

# create DGEList data class
# columns 4 to final column are counts per sample
# columns 1 to 3 are transcript ID, gene name and gene type

y <- DGEList(counts=dat[,4:ncol(dat)], genes=dat[,1:3])

# update the samples table to have patientid and tumourtype columns

y$samples$patientid <- sub("_.*","",rownames(y$samples))
y$samples$tumourtype <- sub(".*_","",rownames(y$samples))
y$samples$tumourtype <- ifelse(y$samples$tumourtype == "Primary","P",ifelse(y$samples$tumourtype == "Recurrent","R",y$samples$tumourtype))

# remove columns corresponding to samples that should be excluded from the analyses

# check all samples occur in metadata

table(y$samples$patientid %in% metadata$Patient.ID)


# filter DGEList using the samples object
# expect to go from 176 samples to 120 samples

keep <- !y$sample$patientid %in% patients.remove # indexes for patient samples not in the patient remove file
y <- y[,keep] # subsetting y


# filter DGEList using the samples object
# expecting to go from 120 samples to 86 samples

keep <- !y$samples$patientid %in% below30
y <- y[,keep]

# TMM normalisation

y <- calcNormFactors(y)

# covert to normalised counts

ynorm <- cpm(y)
rownames(ynorm) <- y$genes$EnsID

# put all expression values in a vector

ylist <- as.vector(ynorm, mode = "numeric")


# remove normalised expression values of zero

ylist <- ylist[ylist > 0.1]

lowerq <- summary(ylist)[2]

lowerq

# now comes the tricky bit
# I am going to split the ynorm object into primary and recurrent samples first

# check that the columns in ynorm are in the same sample order as the rows of y$samples

table(colnames(ynorm) == rownames(y$samples))

idx.primary <- y$samples$tumourtype == "P"
idx.recurrent <- y$samples$tumourtype == "R"

# pull out ynorm columns corresponding to primary tumours

ynormprimary <- ynorm[,idx.primary]

dim(ynormprimary)

# pull out ynorm columns corresponding to recurrent tumours

ynormrecurrent <- ynorm[,idx.recurrent]

dim(ynormrecurrent)

# create new data frame to store % expression at least lower quartile for recurrent and primary tumours for each transcript


func_temp <- function(x){
  x >= lowerq
}

func_temp2 <- function(x){
  ifelse(sum(x)/43 >= 0.2, 1, 0)
}

primtemp <- apply(ynormprimary, 1:2, func_temp)
primtemp2 <- apply(primtemp, 1, func_temp2)

rectemp <- apply(ynormrecurrent, 1:2, func_temp)
rectemp2 <- apply(rectemp, 1, func_temp2)

table(names(primtemp2) == names(rectemp2))

exprres <- cbind(primtemp2, rectemp2)

exprres <- as.data.frame(exprres)
colnames(exprres) <- c("Primary","Recurrent")
exprres$Overall <- ifelse(exprres$Primary == 0 & exprres$Recurrent == 0, "Omit", "Include")

# print list of transcripts to omit

omitidx <- exprres$Overall == "Omit"

length(omitidx) == nrow(y$genes)

omit <- y$genes[omitidx,1]

write.table(omit, "/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/filter3omit.csv")


