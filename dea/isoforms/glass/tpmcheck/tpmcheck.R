library(tidyverse)

# glass #####

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

# plot bar for primary and recurrent tumour from each patient for transcript of interest

transcript.matrix <- as.matrix(transcript.df[,1:2])
transcript.matrix.t <- t(transcript.matrix)
transcript.matrix.t



barplot(transcript.matrix.t, beside = T,
        main = paste0("GLASS TPM Data, Transcript: ",transcript.int),
        xlab = "Patient",
        col = c("red", "green"),
        las = 2)

legend("topright",
       c("Primary", "Recurrent"),
       fill = c("red", "green"))


# in-house #####

tpm.full <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformfpkm_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

tpm <- tpm.full
head(tpm)

patientskeep <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/up-responders/patientskeep.txt", header = F, sep = ",")
head(patientskeep)

# reformat the column names of the dataset for filtering

colnames(tpm) <- sub("Primary","P",colnames(tpm))
colnames(tpm) <- sub("Recurrent","R",colnames(tpm))
colnames(tpm) <- sub("_FPKM","",colnames(tpm))
colnames(tpm)

keep <- sub(".{2}$","",colnames(tpm)) %in% patientskeep[,1]
keep[1] <- TRUE

tpm <- tpm[,keep]

# get list of top DE isoforms in up-responders from Stead data

stead.up <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv",
                       header = T, sep = ",")
stead.up$EnsID <- sub("\\..*","",stead.up$Row.names)

# list of patients

patients <- unique(sub(".{2}$","",colnames(tpm)[2:ncol(tpm)]))

# transcript of interest

transcript.int <- "ENST00000335585"

transcript <- tpm[sub(".{2}$","",tpm$EnsID) == transcript.int,]

transcript.df <- as.data.frame(matrix(ncol=2, nrow=(ncol(transcript)-1)/2))
rownames(transcript.df) <- patients
colnames(transcript.df) <- c("Primary", "Recurrent")

for (i in 1:nrow(transcript.df)){
  idx <- match(rownames(transcript.df)[i],sub("_P$","",colnames(transcript)))
  transcript.df[i,1] <- transcript[1,idx]
  idx <- match(rownames(transcript.df)[i],sub("_R$","",colnames(transcript)))
  transcript.df[i,2] <- transcript[1,idx]
}

# plot bar for primary and recurrent tumour from each patient for transcript of interest

transcript.matrix <- as.matrix(transcript.df[,1:2])
transcript.matrix.t <- t(transcript.matrix)
transcript.matrix.t

barplot(transcript.matrix.t, beside = T,
        main = paste0("Stead FPKM Data, Transcript: ",transcript.int),
        col = c("darkblue", "lightblue"),
        las = 2)

legend("topright",
       c("Primary", "Recurrent"),
       fill = c("darkblue", "lightblue"))









    