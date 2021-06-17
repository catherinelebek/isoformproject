# import fpkm data

# datfull.fpkm <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/PvR_isoformfpkm_all.txt", header = TRUE)
datfull.fpkm <- read.delim("/nobackup/bs20chlb/inputdata/PvR_isoformfpkm_all.txt", header = TRUE)

dat.fpkm <- datfull.fpkm

# import metadata

# metadata <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/Metadata.csv", header = TRUE)
metadata <- read.csv("/nobackup/bs20chlb/inputdata/Metadata.csv", header = TRUE)

# import list of patients to remove based on metadata values

# patients.remove <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/patients_remove.txt", header = FALSE)
patients.remove <- read.delim("/nobackup/bs20chlb/inputdata/patients_remove.txt", header = FALSE)

# convert to vector

patients.remove <- as.vector(t(patients.remove))

# import list of patients to remove based on reads < 30m

# below30 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/below30.txt", header = FALSE)
below30 <- read.delim("/nobackup/bs20chlb/inputdata/below30.txt", header = FALSE)

# covert datamframe to vector

below30 <- as.vector(t(below30))

# rearrange columns

c0 <- grep("EnsID", colnames(dat.fpkm))
c1 <- grep("GeneName", colnames(dat.fpkm))
c2 <- grep("GeneType", colnames(dat.fpkm))
c3 <- c1 - 1
c4 <- c2 + 1
rownames(dat.fpkm) <- dat.fpkm[,c0]
dat.fpkm <- dat.fpkm[,c(2:c3,c4:ncol(dat.fpkm))]
colnames(dat.fpkm)

# remove samples based on metadata values

keep <- !sub("_.*","",colnames(dat.fpkm[1:ncol(dat.fpkm)])) %in% patients.remove
dat.fpkm <- dat.fpkm[,keep]

# remove samples based on read counts <3m

keep <- !sub("_.*","",colnames(dat.fpkm[1:ncol(dat.fpkm)])) %in% below30
dat.fpkm <- dat.fpkm[,keep]

# remove fpkm from column names

colnames(dat.fpkm) <- sub("_FPKM","",colnames(dat.fpkm))
colnames(dat.fpkm) <- sub("Primary","P",colnames(dat.fpkm))
colnames(dat.fpkm) <- sub("Recurrent","R",colnames(dat.fpkm))


# check all samples occur in metadata

table(sub("_.*","",colnames(dat.fpkm)) %in% metadata$Patient.ID)

# put all expression values in a vector

ylist <- as.numeric(unlist(dat.fpkm))
ylist <- as.vector(ylist, mode = "numeric")

# remove normalised expression values of zero

ylist <- ylist[ylist != 0]

# pull out the lower quartile

lowerq <- summary(ylist)[2]

lowerq

# now comes the tricky bit
# I am going to split the dat.fpkm object into primary and recurrent samples first

idx.primary <- "P" == sub(".*_","",colnames(dat.fpkm))
idx.recurrent <- "R" == sub(".*_","",colnames(dat.fpkm))

# pull out ynorm columns corresponding to primary tumours

datfpkmprimary <- dat.fpkm[,idx.primary]

# pull out ynorm columns corresponding to recurrent tumours

datfpkmrecurrent <- dat.fpkm[,idx.recurrent]

# create new data frame to store % expression at least lower quartile for recurrent and primary tumours for each transcript

exprres <- data.frame()

for (i in 1:nrow(datfpkmprimary)){ # for every transcript
  temp <- c()
  for (j in 1:ncol(datfpkmprimary)){ # take every patient
    temp[j] <- datfpkmprimary[i,j] >= lowerq # determine if the expression is higher or equal to the lower quartile
  }
  res <- sum(temp) / ncol(datfpkmprimary) # determine what % of all patiets have expression higher than or equal to the lower quartile
  exprres[i,1] <- ifelse(res >= 0.2, 1, 0) # determine is this % is great than or equal to 20%
}


for (i in 1:nrow(datfpkmrecurrent)){ # for every transcript
  temp <- c()
  for (j in 1:ncol(datfpkmrecurrent)){ # take every patient
    temp[j] <- datfpkmrecurrent[i,j] >= lowerq # determine if the expression is higher or equal to the lower quartile
  }
  res <- sum(temp) / ncol(datfpkmrecurrent) # determine what % of all patiets have expression higher than or equal to the lower quartile
  exprres[i,2] <- ifelse(res >= 0.2, 1, 0) # determine is this % is great than or equal to 20%
}

colnames(exprres) <- c("Primary","Recurrent")
exprres$Overall <- ifelse(exprres$Primary == 0 & exprres$Recurrent == 0, "Omit", "Include")

table(exprres$Overall)


# print list of transcripts to omit

omitidx <- exprres$Overall == "Omit"

length(omitidx) <- y$genes[omitidx,1]

omit <- rownames(dat.fpkm)[omitidx]

# write.table(omit, "/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/lowexpressionomit.csv")

write.table(omit, "/nobackup/bs20chlb/outputdata/fpkmlowexpressionomit.csv")
