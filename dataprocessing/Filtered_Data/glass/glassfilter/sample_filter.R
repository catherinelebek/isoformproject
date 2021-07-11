glassdatfull <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/glassdata/transcript_count_matrix_all_samples.tsv", header = T, sep = "\t")

glassdat <- glassdatfull

patientskeep <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/glassdata/glasspatientkeep.csv", header = T, sep = ",")

# reformat the column names of the glass dataset for filtering

colnames(glassdat) <- sub("...R.RNA.*","",colnames(glassdat))
colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)
colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)
colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)

colnames(glassdat)

patientskeep <- c(patientskeep[,2], patientskeep[,3])

keep <- colnames(glassdat) %in% patientskeep
keep[1] <- TRUE
keep[2] <- TRUE

glassdat <- glassdat[,keep]

table(colnames(glassdat))

head(glassdat)

# import list of transcripts to remove

lowexp <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/isoforms/glass/glassfilter/transcriptsomit.csv",
                     header = T, sep = ",")

lowexp <- c(t(lowexp))

head(glassdat)
# filter out the transcripts with low expression values

keep <- !glassdat[,1] %in% lowexp
glassdat <- glassdat[keep,] # this is the final set of data


write.table(glassdat, "/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/glass/glassfilter/PvR_isoformCounts_filtered.txt",
            sep = "\t")
