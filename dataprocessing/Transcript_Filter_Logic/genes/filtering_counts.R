# this script gives a list of genes we want to keep based on count data from patients we want to include

# import libraries

library(edgeR)

datfull.counts <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/PvR_geneCounts_all_LS_23062021.txt.txt", 
                             header = TRUE, sep = "\t")

dat <- datfull.counts

# import list of patients to remove based on metadata values

patientskeep <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/patientskeep.txt", header = FALSE)

# convert to vector

patientskeep <- as.vector(t(patientskeep))

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

y$samples$Patient.ID <- sub("rimary","",rownames(y$samples))
y$samples$Patient.ID <- sub("ecurrent","",y$samples$Patient.ID)
y$samples$Patient.ID <- gsub(".{2}$","",y$samples$Patient.ID)

y$samples$tumour.type <- sub("rimary","",rownames(y$samples))
y$samples$tumour.type <- sub("ecurrent","",y$samples$tumour.type)
y$samples$tumour.type <- gsub(".*_","",y$samples$tumour.type)



# remove columns corresponding to samples that should be excluded from the analyses

# filter DGEList using the samples object

# expect to go from 240 samples to 66x2 (132) samples

keep <- y$samples$Patient.ID %in% patientskeep # indexes for patient samples in the patientkeep file
y <- y[,keep] # subsetting y

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

idx.primary <- y$samples$tumour.type == "P"
idx.recurrent <- y$samples$tumour.type == "R"

# pull out ynorm columns corresponding to primary tumours

ynormprimary <- ynorm[,idx.primary]

dim(ynormprimary)

# pull out ynorm columns corresponding to recurrent tumours

ynormrecurrent <- ynorm[,idx.recurrent]

dim(ynormrecurrent)

# create new data frame to store % expression at least lower quartile for recurrent and primary tumours for each transcript

n <- 66

func_temp <- function(x){
  x >= lowerq
}

func_temp2 <- function(x){
  ifelse(sum(x)/n >= 0.2, 1, 0)
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

write.table(omit, "/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/genes/transcriptsomit.csv", row.names = F)


