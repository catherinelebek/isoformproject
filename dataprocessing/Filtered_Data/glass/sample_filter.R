patientskeep <- read.delim("~/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/isoforms/glassdata/samplefilters/glasspatientkeep.csv", header = T, sep = ",")

glassdatfull <- read.delim("~/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/isoforms/glassdata/transcript_count_matrix_all_samples.tsv", header = T, sep = "\t")

head(glassdat)
colnames(glassdat)

glassdat <- glassdatfull

# reformat the column names of the glass dataset for filtering

colnames(glassdat) <- sub("...R.RNA.*","",colnames(glassdat))


colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)
colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)
colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)

colnames(glassdat)

patientskeep <- c(patientskeep[,2], patientskeep[,3])

keep <- colnames(glassdat) %in% patientskeep
keep <- c(TRUE,TRUE,keep)

glassdat <- glassdat[,keep]

table(colnames(glassdat))

write.table(glassdat, "/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/isoforms/glass/PvR_isoformCounts_filtered.txt",
            sep = "\t")
