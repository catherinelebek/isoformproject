# import libraries

library(edgeR)

# import count data

datfull <- read.delim("/nobackup/bs20chlb/inputdata/PvR_isoformCounts_all.txt", header = TRUE)
dat <- datfull

# import metadata

metadata <- read.csv("/nobackup/bs20chlb/inputdata/Metadata.csv", header = TRUE)

# rearrange columns

c0 <- grep("EnsID", colnames(dat))
c1 <- grep("GeneName", colnames(dat))
c2 <- grep("GeneType", colnames(dat))
c3 <- c1 - 1
c4 <- c2 + 1
dat <- dat[,c(c0,c1,c2,2:c3,c4:ncol(dat))]
colnames(dat)

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

# import list of patients to remove

patients.remove <- read.delim("/nobackup/bs20chlb/inputdata/patients_remove.txt", header = FALSE)

# convert list to vector

patients.remove <- as.vector(t(patients.remove))

# filter DGEList using the samples object

keep <- !y$sample$patientid %in% patients.remove # indexes for patient samples not in the patient remove file
y <- y[,keep] # subsetting y

# TMM normalisation

y <- calcNormFactors(y)

# covert to normalised counts

ynorm <- cpm(y)

# put all expression values in a vector

ylist <- as.vector(ynorm, mode = "numeric")

# remove normalised expression values of zero

ylist <- ylist[ylist != 0]

# pull out the lower quartile

lowerq <- summary(ylist)[2]

# now comes the tricky bit
# I am going to split the ynorm object into primary and recurrent samples first

# check that the columns in ynorm are in the same sample order as the rows of y$samples

table(colnames(ynorm) == rownames(y$samples))

idx.primary <- y$samples$tumourtype == "P"
idx.recurrent <- y$samples$tumourtype == "R"

# pull out ynorm columns corresponding to primary tumours

ynormprimary <- ynorm[,idx.primary]

# pull out ynorm columns corresponding to recurrent tumours

ynormrecurrent <- ynorm[,idx.recurrent]

# create new data frame to store % expression at least lower quartile for recurrent and primary tumours for each transcript

exprres <- data.frame()

for (i in 1:nrow(ynormprimary)){
  temp <- c()
  for (j in 1:ncol(ynormprimary)){
    temp[j] <- ynormprimary[i,j] >= lowerq
  }
  temp <- sum(temp) / ncol(ynormprimary)
  exprres[i,1] <- ifelse(temp >= 0.2, 1, 0)
}


for (i in 1:nrow(ynormrecurrent)){ # for every transcript
  temp <- c()
  for (j in 1:ncol(ynormrecurrent)){ # take every patient
    temp[j] <- ynormrecurrent[i,j] >= lowerq # determine if the expression is higher or equal to the lower quartile
  }
  temp <- sum(temp) / ncol(ynormrecurrent) # determine what % of all patnets have expression higher than or equal to the lower quartile
  exprres[i,2] <- ifelse(temp >= 0.2, 1, 0)
}

colnames(exprres) <- c("Primary","Recurrent")
exprres$Overall <- ifelse(exprres$Primary == 0 & exprres$Recurrent == 0, "Omit", "Include")
exprres

keep <- exprres$Overall == "Include"

# check that keep is the same length as y

length(keep) == nrow(y$counts)

# print list of transcripts to omit
exprres
omitidx <- exprres$Overall == "Omit"
omit <- y$genes[omitidx,1]

write.table(omit, "/nobackup/bs20chlb/outputdata/lowexpressionomit.csv")


