library(stringr)

# import count data

datfull <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/PvR_geneCounts_all_LS_23062021.txt.txt",
                      header = TRUE, sep = "\t")


dat <- datfull


# import list of patients to remove based on metadata values

patientskeep <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/up-responders/patientskeep.txt", header = FALSE)
# patients.remove <- read.delim("/nobackup/bs20chlb/inputdata/patients_remove.txt", header = FALSE)

# convert to vector

patientskeep <- as.vector(t(patientskeep))


# rearrange columns

c0 <- grep("EnsID", colnames(dat))
c1 <- grep("GeneName", colnames(dat))
c2 <- grep("GeneType", colnames(dat))
c3 <- c0 + 1
c4 <- c1 - 1
c5 <- c2 + 1
dat <- dat[,c(c0,c1,c2,c3:c4,c5:ncol(dat))]
colnames(dat)

# ensure P and R

colnames(dat) <- sub("Primary","P",colnames(dat))
colnames(dat) <- sub("Recurrent","R",colnames(dat))

# remove samples based on metadata values
# do this by creating a vector of column names in the same order as in dat

cols <- colnames(dat[,4:ncol(dat)])
cols <- str_sub(cols, start = 1, end = -3)
cols 

keep <- cols %in% patientskeep
add <- c(TRUE,TRUE,TRUE)
keep <- c(add,keep)

dat <- dat[,keep]

head(dat)

write.table(dat, "/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/genes/up-responders/PvR_geneCounts_filtered.txt",
            sep = "\t")


