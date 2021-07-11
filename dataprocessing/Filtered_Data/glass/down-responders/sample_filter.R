glassdatfull <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/glassdata/transcript_count_matrix_all_samples.tsv", header = T, sep = "\t")

glassdat <- glassdatfull

patientskeep <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/glassdata/glasspatientkeep.csv", header = T, sep = ",")

jarid <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/glassdata/JARID2_results.txt", header = T, sep = "\t")

# get responder type status

jarid$NES <- as.numeric(jarid$NES)
jarid$direction <- ifelse(is.na(jarid$NES), jarid$ES, jarid$NES)
jarid$responder.type <- ifelse(jarid$direction > 0, "Up", "Down")
jarid

# reformat jarid so . is -

jarid$Patient <- sub(".","-",jarid$Patient,fixed=TRUE)
jarid$Patient <- sub(".","-",jarid$Patient,fixed=TRUE)

# reformat the column names of the glass dataset for filtering

colnames(glassdat) <- sub("...R.RNA.*","",colnames(glassdat))
colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)
colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)
colnames(glassdat) <- sub(".","-",colnames(glassdat), fixed = TRUE)

colnames(glassdat)

# pull down-responders from patients keep - 9 patients

down <- patientskeep[,1] %in% jarid[jarid$responder.type == "Down", 1]
patientskeep <- patientskeep[down,]

patientskeep <- c(patientskeep[,2], patientskeep[,3])

keep <- colnames(glassdat) %in% patientskeep
keep[1] <- TRUE
keep[2] <- TRUE

glassdat <- glassdat[,keep]

table(colnames(glassdat))

# import list of isoforms analysed in up-responders from the stead data

stead.down <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results.csv", header = T, sep = ",")
stead.down <- stead.down$Row.names

stead.down <- sub("\\..*","",stead.down)

# filter out the transcripts with low expression values

keep <- glassdat[,1] %in% stead.down
glassdat <- glassdat[keep,] # this is the final set of data



write.table(glassdat, "/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/glass/down-responders/PvR_isoformCounts_filtered.txt",
            sep = "\t")
