# this script gives a list of transcripts we want to keep based on count data from patients we want to include

# import libraries

library(edgeR)

datfull.counts <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/glassdata/transcript_count_matrix_all_samples.tsv", 
                             header = TRUE, sep = "\t")

dat <- datfull.counts

head(dat)

# import list of patients to remove based on metadata values

patientskeep <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/glassdata/glasspatientkeep.csv", 
                           header = TRUE, sep = ",")

# create DGEList data class
# columns 4 to final column are counts per sample
# columns 1 to 3 are transcript ID, gene name and gene type

y <- DGEList(counts=dat[,3:ncol(dat)], genes=dat[,1:2])

# update the samples table to have patientid and tumourtype columns

y$samples$Patient.ID <- sub("...R.RNA.*","",rownames(y$samples))
y$samples$Patient.ID <- sub(".","-",y$samples$Patient.ID, fixed = TRUE)
y$samples$Patient.ID  <- sub(".","-",y$samples$Patient.ID, fixed = TRUE)
y$samples$Patient.ID <- sub(".","-",y$samples$Patient.ID, fixed = TRUE)
y$samples$tumour.type <- sub(".*-","",y$samples$Patient.ID)


head(y$samples)


# remove columns corresponding to samples that should be excluded from the analyses

# filter DGEList using the samples object

patientskeep <- c(t(patientskeep[,2]),t(patientskeep[,3]))
keep <- y$samples$Patient.ID %in% patientskeep # indexes for patient samples in the patientkeep file
y <- y[,keep] # subsetting y

# TMM normalisation

y <- calcNormFactors(y)

# covert to normalised counts

ynorm <- cpm(y)
rownames(ynorm) <- y$genes$target_id

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

idx.primary <- y$samples$tumour.type == "R1"
idx.recurrent <- y$samples$tumour.type == "TP"

# pull out ynorm columns corresponding to primary tumours

ynormprimary <- ynorm[,idx.primary]

dim(ynormprimary)

# pull out ynorm columns corresponding to recurrent tumours

ynormrecurrent <- ynorm[,idx.recurrent]

dim(ynormrecurrent)

# create new data frame to store % expression at least lower quartile for recurrent and primary tumours for each transcript

n <- ncol(ynormprimary)

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

write.table(omit, "/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/glass/glassfilter/transcriptsomit.csv", row.names = F)


