# import fpkm data

datfull <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/PvR_isoformfpkm_all.txt", header = TRUE)
# datfull.fpkm <- read.delim("/nobackup/bs20chlb/inputdata/PvR_isoformfpkm_all.txt", header = TRUE)

dat <- datfull

# import metadata

metadata <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/samplefilters/Metadata.csv", header = TRUE)
# metadata <- read.csv("/nobackup/bs20chlb/inputdata/Metadata.csv", header = TRUE)

# import list of patients to remove based on metadata values

patients.remove <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/samplefilters/patients_remove.txt", header = FALSE)
# patients.remove <- read.delim("/nobackup/bs20chlb/inputdata/patients_remove.txt", header = FALSE)

# convert to vector

patients.remove <- as.vector(t(patients.remove))

# import list of patients to remove based on reads < 30m

below30 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/samplefilters/below30.txt", header = FALSE)
# below30 <- read.delim("/nobackup/bs20chlb/inputdata/below30.txt", header = FALSE)

# covert datamframe to vector

below30 <- as.vector(t(below30))

# rearrange columns

c0 <- grep("EnsID", colnames(dat))
c1 <- grep("GeneName", colnames(dat))
c2 <- grep("GeneType", colnames(dat))
c3 <- c0 + 1
c4 <- c1 - 1
c5 <- c2 + 1
dat <- dat[,c(c0,c1,c2,c3:c4,c5:ncol(dat))]
colnames(dat)

# remove samples based on metadata values

keep <- !sub("_.*","",colnames(dat)) %in% patients.remove
dat <- dat[,keep]

# remove samples based on read counts <3m

keep <- !sub("_.*","",colnames(dat)) %in% below30
dat <- dat[,keep]

# ensure P and R

colnames(dat) <- sub("Primary","P",colnames(dat))
colnames(dat) <- sub("Recurrent","R",colnames(dat))

# check all samples occur in metadata

table(sub("_.*","",colnames(dat[4:ncol(dat)])) %in% metadata$Patient.ID)

lowexp <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/nonstringentinclude.txt", header = TRUE, sep = "\t")
# lowexp <- read.csv("/nobackup/bs20chlb/inputdata/lowexpressionomit.csv", header = TRUE)

# convert from data frame to vector
lowexp <- as.vector(t(lowexp[,1]))

# filter out the transcripts with low expression values

keep <- dat[,1] %in% lowexp
dat <- dat[keep,] # this is the final set of data

# remove fpkm from column names

colnames(dat) <- sub("_FPKM","",colnames(dat))
head(dat)

# write to file

write.table(dat, "/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/nonstringent/PvR_isoformfpkm_filtered.txt")

