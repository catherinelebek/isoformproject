# this script gives a list of transcripts we want to keep based on count data from patients we want to include

# import libraries

library(edgeR)
library(tidyverse)

datfull.counts <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt", 
                             header = TRUE, sep = "\t")

dat <- datfull.counts
colnames(dat) <- sub("rimary","",colnames(dat))
colnames(dat) <- sub("ecurrent","", colnames(dat))

# import list of patients to remove based on metadata values

patientskeep <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/patientskeep_batch.txt", header = FALSE)

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

y$samples$Patient.ID <- gsub(".{2}$","",rownames(y$samples))
y$samples$tumour.type <- gsub(".*_","",rownames(y$samples))

# update the samples table to include whether an up or down responder

metadata <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/MetaData_GT_250621.txt", 
                      header = TRUE, sep = "\t")

merge <- metadata %>% select(Patient.ID, NES)

y$samples <- left_join(y$samples, merge, by = c("Patient.ID" = "Patient.ID"))


# remove columns corresponding to samples that should be excluded from the analyses

# filter DGEList using the samples object

keep <- paste(y$samples$Patient.ID,"_",y$samples$tumour.type,sep="") %in% patientskeep # indexes for patient samples in the patientkeep file
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

table(colnames(ynorm) == paste(y$samples$Patient.ID,"_",y$samples$tumour.type,sep=""))

idx.up <- y$samples$NES > 0
idx.down <- y$samples$NES < 0

# pull out ynorm columns corresponding to primary tumours

ynormup<- ynorm[,idx.up]

dim(ynormup)

# pull out ynorm columns corresponding to recurrent tumours

ynormdown <- ynorm[,idx.down]

dim(ynormdown)

# create new data frame to store % expression at least lower quartile for recurrent and primary tumours for each transcript

n.up <- ncol(ynormup)
n.down <- ncol(ynormdown)

func_temp <- function(x){
  x >= lowerq
}

func_temp2 <- function(x){
  ifelse(sum(x)/n.up >= 0.2, 1, 0)
}

func_temp3 <- function(x){
  ifelse(sum(x)/n.down >= 0.2, 1, 0)
}


uptemp <- apply(ynormup, 1:2, func_temp)
uptemp2 <- apply(uptemp, 1, func_temp2)

downtemp <- apply(ynormdown, 1:2, func_temp)
downtemp2 <- apply(downtemp, 1, func_temp3)

table(names(uptemp2) == names(downtemp2))

exprres <- cbind(uptemp2, downtemp2)

exprres <- as.data.frame(exprres)
colnames(exprres) <- c("Up","Down")
exprres$Overall <- ifelse(exprres$Up == 0 & exprres$Down == 0, "Omit", "Include")

# print list of transcripts to omit

omitidx <- exprres$Overall == "Omit"

length(omitidx) == nrow(y$genes)

omit <- y$genes[omitidx,1]

write.table(omit, "/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/transcriptsomit_batch.csv", row.names = F)


