library(tidyverse)

tpm.full <- readr::read_tsv("~/Documents/Semester3/Project/InputData/isoforms/glassdata/transcript_tpm_matrix_all_samples.tsv.gz")

tpm <- tpm.full

patientskeep <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/glassdata/glasspatientkeep.csv", header = T, sep = ",")

jarid <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/glassdata/JARID2_results.txt", header = T, sep = "\t")

# get responder type status

jarid$NES <- as.numeric(jarid$NES)
jarid$direction <- ifelse(is.na(jarid$NES), jarid$ES, jarid$NES)
jarid$responder.type <- ifelse(jarid$direction > 0, "Up", "Down")
jarid
table(jarid$responder.type)

# reformat jarid so . is -

jarid$Patient <- sub(".","-",jarid$Patient,fixed=TRUE)
jarid$Patient <- sub(".","-",jarid$Patient,fixed=TRUE)

# reformat the column names of the glass dataset for filtering

colnames(tpm) <- sub("-..R-RNA-.*","",colnames(tpm))

# pull up-responders from patients keep - 14 patients

up <- patientskeep[,1] %in% jarid[jarid$responder.type == "Up", 1]
patientskeep <- patientskeep[up,]

patientskeep <- c(patientskeep[,2], patientskeep[,3])

keep <- colnames(tpm) %in% patientskeep
keep[1] <- TRUE
keep[2] <- TRUE

tpm <- tpm[,keep]

# get list of top DE isoforms in up-responders from Stead data

stead.up <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv",
                       header = T, sep = ",")
stead.up$EnsID <- sub("\\..*","",stead.up$Row.names)

# list of patients

patients <- unique(sub(".{3}$","",colnames(tpm)[3:ncol(tpm)]))

# transcript of interest

transcript.int <- "ENST00000335585"

transcript <- tpm[tpm$target_id == transcript.int,]
transcript <- as.data.frame(transcript)
head(transcript)

transcript.df <- as.data.frame(matrix(ncol=2, nrow=(ncol(transcript)-2)/2))
rownames(transcript.df) <- patients
colnames(transcript.df) <- c("Primary", "Recurrent")

for (i in 1:nrow(transcript.df)){
    idx <- match(rownames(transcript.df)[i],sub("-TP","",colnames(transcript)))
    transcript.df[i,1] <- transcript[1,idx]
    idx <- match(rownames(transcript.df)[i],sub("-R1","",colnames(transcript)))
    transcript.df[i,2] <- transcript[1,idx]
}



transcript.df$LFC <- log2(transcript.df$Recurrent) - log2(transcript.df$Primary)
boxplot(transcript.df$LFC)
transcript.df 

par(mfrow = c(1,2))
boxplot(transcript.df$Primary)
boxplot(transcript.df$Recurrent)

head(transcript)
match(rownames(transcript.df)[1],sub("-TP","",colnames(tpm)))

    