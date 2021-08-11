library(tidyverse)
library(reshape2)

patienttype <- "down-responders/"
type <- "Down"

# glass #####

tpm.full <- readr::read_tsv("~/Documents/Semester3/Project/InputData/isoforms/glassdata/transcript_tpm_matrix_all_samples.tsv.gz")

tpm.glass <- tpm.full

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

colnames(tpm.glass) <- sub("-..R-RNA-.*","",colnames(tpm.glass))

# pull up-responders from patients keep - 14 patients

up <- patientskeep[,1] %in% jarid[jarid$responder.type == type, 1]
patientskeep <- patientskeep[up,]

patientskeep <- c(patientskeep[,2], patientskeep[,3])

keep <- colnames(tpm.glass) %in% patientskeep
keep[1] <- TRUE
keep[2] <- TRUE

tpm.glass <- tpm.glass[,keep]

# list of patients

patients <- unique(sub(".{3}$","",colnames(tpm.glass)[3:ncol(tpm.glass)]))

# transcript of interest

transcript.int <- "ENST00000315927"

transcript <- tpm.glass[tpm.glass$target_id == transcript.int,]
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

transcript.matrix.glass <- as.matrix(transcript.df[,1:2])

transcript.matrix.glass <- melt(transcript.matrix.glass[,c("Primary","Recurrent")],id.vars = 1)

p <- ggplot(data = transcript.matrix.glass, aes(x=Var1, y=value, fill=Var2)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Isoform Count") +
  theme(axis.text.x = element_text(size = "12", angle = 90), 
        axis.text.y = element_text(size = "12"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  ggpubr::rremove("grid")

p

# in-house #####

tpm.full.inhouse <- read.delim("~/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformfpkm_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

tpm <- tpm.full.inhouse
head(tpm)

patientskeep <- read.delim(paste0("~/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/",patienttype,"patientskeep.txt"), header = F, sep = ",")
head(patientskeep)

# reformat the column names of the dataset for filtering

colnames(tpm) <- sub("Primary","P",colnames(tpm))
colnames(tpm) <- sub("Recurrent","R",colnames(tpm))
colnames(tpm) <- sub("_FPKM","",colnames(tpm))
colnames(tpm)

keep <- sub(".{2}$","",colnames(tpm)) %in% patientskeep[,1]
keep[1] <- TRUE

tpm <- tpm[,keep]

# list of patients

patients <- unique(sub(".{2}$","",colnames(tpm)[2:ncol(tpm)]))

# transcript of interest

transcript.int <- "ENST00000377819"

transcript <- tpm[sub("\\..*","",tpm$EnsID) == transcript.int,]

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

transcript.matrix <- melt(transcript.matrix[,c("Primary","Recurrent")],id.vars = 1)

p <- ggplot(data = transcript.matrix, aes(x=Var1, y=value, fill=Var2)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Isoform Count") +
  theme(axis.text.x = element_text(size = "12", angle = 90), 
        axis.text.y = element_text(size = "12"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  ggpubr::rremove("grid")

p


transcript.matrix

# merge data

transcript.matrix.master <- rbind(transcript.matrix.glass, transcript.matrix)
transcript.matrix.master$source <- ifelse(gsub("-.*","",transcript.matrix.master$Var1) == "GLSS"
                                          | gsub("-.*","",transcript.matrix.master$Var1) == "TCGA",
                                          "GLASS","In-house")

transcript.matrix.master$legend <- paste(transcript.matrix.master$Var2, transcript.matrix.master$source)
table(transcript.matrix.master$source)

p <- ggplot(data = transcript.matrix.master, aes(x=Var1, y=value, fill= legend)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Isoform Count") +
  theme(axis.text.x = element_text(size = "12", angle = 90), 
        axis.text.y = element_text(size = "12"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  ggpubr::rremove("grid")

p


# running analysis on UvD responders #####

tpm.glass <- tpm.full

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

colnames(tpm.glass) <- sub("-..R-RNA-.*","",colnames(tpm.glass))

# pull up-responders from patients keep - 14 patients

up <- patientskeep[,1] %in% jarid[jarid$responder.type == "Down", 1]
patientskeep <- patientskeep[up,]

# pull only primary tumour expression

patientskeep <- patientskeep[,2]

keep <- colnames(tpm.glass) %in% patientskeep
keep[1] <- TRUE
keep[2] <- TRUE

tpm.glass <- tpm.glass[,keep]

# list of patients

patients <- unique(sub(".{3}$","",colnames(tpm.glass)[3:ncol(tpm.glass)]))

# transcript of interest

transcript.int <- "ENST00000554738"

transcript <- tpm.glass[tpm.glass$target_id == transcript.int,]
transcript <- as.data.frame(transcript)
head(transcript)

transcript.df <- as.data.frame(matrix(ncol=2, nrow=(ncol(transcript)-2)))
transcript.df[,1] <- colnames(transcript)[c(-1,-2)]
colnames(transcript.df) <- c("Patient","Up")

for (i in 1:nrow(transcript.df)){
  idx <- match(transcript.df[i,1],colnames(transcript))
  transcript.df[i,2] <- transcript[1,idx]
}

transcript.df

# plot bar 

p <- ggplot(data = transcript.df, aes(x=Patient, y=Up)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Isoform Count") +
  theme(axis.text.x = element_text(size = "12", angle = 90), 
        axis.text.y = element_text(size = "12"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  ggpubr::rremove("grid")

p 

p1 <- ggplot(data = transcript.df, aes(x=Patient, y=Up)) +
  geom_col(position = "dodge") +
  theme_bw() +
  scale_x_discrete(name = "Patient") +
  scale_y_continuous(name = "Isoform Count") +
  theme(axis.text.x = element_text(size = "12", angle = 90), 
        axis.text.y = element_text(size = "12"),
        axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
        axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
        plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
  ggpubr::rremove("grid")

p1

transcript









    