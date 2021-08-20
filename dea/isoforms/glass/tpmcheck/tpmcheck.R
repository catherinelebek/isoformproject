library(tidyverse)
library(reshape2)
library(patchwork)

patienttype <- "up-responders/"
type <- "Up"

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
transcript.int <- "ENST00000545061" # SCN8A-203
transcript.int <- "ENST00000354534" # SCN8A-201
transcript.int <- "ENST00000377819" # RTN3-206
transcript.int <- "ENST00000626839" # KCNQ2-213

transcript <- tpm.glass[tpm.glass$target_id == transcript.int,]
transcript <- as.data.frame(transcript)
head(transcript)

transcript.df <- as.data.frame(matrix(ncol=3, nrow=(ncol(transcript)-2)/2))
rownames(transcript.df) <- patients
colnames(transcript.df) <- c("Primary", "Recurrent","Patient.Num")

for (i in 1:nrow(transcript.df)){
    idx <- match(rownames(transcript.df)[i],sub("-TP","",colnames(transcript)))
    transcript.df[i,1] <- transcript[1,idx]
    idx <- match(rownames(transcript.df)[i],sub("-R1","",colnames(transcript)))
    transcript.df[i,2] <- transcript[1,idx]
    transcript.df[i,3] <- i
}

transcript.df

# plot bar for primary and recurrent tumour from each patient for transcript of interest

transcript.matrix.glass <- transcript.df[,c(3,1,2)]
transcript.matrix.glass <- as.matrix(transcript.matrix.glass[,1:3])

transcript.matrix.glass <- melt(transcript.matrix.glass[,c("Primary","Recurrent")], id.vars = 1)
transcript.matrix.glass$Patient.Num <- seq(1,14,1)

patients <- as.character(seq(1,14,1))

# SCN8A-201, padj = 9.32e-4
# RTN3-206, padj = 0.14

p <- ggplot(data = transcript.matrix.glass, aes(x=Var1, y=value)) +
  geom_col(position = "dodge", aes(fill = Var2, color = Var2)) +
  theme_bw() +
  scale_x_discrete(name = "Patient", labels = patients) +
  scale_y_continuous(name = "TPM") +
  ggtitle("RTN3-206, padj = 0.14") +
  scale_fill_manual(values = c("Primary" = "#190F9733", "Recurrent" = "#208C1733")) +
  scale_color_manual(values = c("Primary" = "#190F97", "Recurrent" = "#208C17")) +
  theme(axis.text.x = element_text(size = "10", angle = 90, hjust = 0.95), 
        axis.text.y = element_text(size = "12"),
        axis.title.x = element_text(size = "12", margin = margin(t=10)),
        axis.title.y = element_text(size = "12", margin = margin(r=5)),
        plot.title = element_text(face = "bold", size=12),
        legend.position = "none") +
  ggpubr::rremove("grid")

p

# create vector to store plots

plots_paired_dea <- vector(mode = "list", length = 4)

plots_paired_dea[[3]] <- p

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

patients <- as.character(seq(1,44,1))

# SCN8A-201, padj = 1.49e-31
# RTN3-206, padj = 2.78e-16

p2 <- ggplot(data = transcript.matrix, aes(x=Var1, y=value)) +
  geom_col(position = "dodge", aes(fill = Var2, color = Var2)) +
  theme_bw() +
  scale_x_discrete(name = "Patient", labels = patients) +
  scale_y_continuous(name = "FPKM") +
  ggtitle("RTN3-206, padj = 2.78e-16") +
  scale_fill_manual(values = c("Primary" = "#190F9733", "Recurrent" = "#208C1733")) +
  scale_color_manual(values = c("Primary" = "#190F97", "Recurrent" = "#208C17")) +
  labs(fill = "Tumour Type", color = "Tumour Type") +
  theme(axis.text.x = element_text(size = "10", angle = 90, hjust = 0.95), 
        axis.text.y = element_text(size = "12"),
        axis.title.x = element_text(size = "12", margin = margin(t=10)),
        axis.title.y = element_text(size = "12", margin = margin(r=5)),
        plot.title = element_text(face = "bold", size=12),
        legend.position = "none") +
  ggpubr::rremove("grid")

p2


plots_paired_dea[[4]] <- p2

(plots_paired_dea[[1]] | plots_paired_dea[[3]]) /
    plots_paired_dea[[2]] /
    plots_paired_dea[[4]]

# running analysis on UvD responders #####

tpm.glass <- tpm.full

# load full list of transcripts in order to pull through gene names
genelist <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")
# extract just transcripts and gene names
genelist <- genelist[,c(1,10)]
genelist$EnsID <- sub("\\..*","",genelist$EnsID)

# add gene names to TPM data

tpm.glass <- left_join(tpm.glass, genelist, by = c("target_id" = "EnsID"))
tpm.glass <- as.data.frame(tpm.glass)

# bring gene name forward as a column

tpm.glass <- tpm.glass[,c(1,387,2,3:386)]

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

# pull one responder type from patients keep

up <- patientskeep[,1] %in% jarid[jarid$responder.type == "Up", 1]
down <- patientskeep[,1] %in% jarid[jarid$responder.type == "Down", 1]

patientskeep.up <- patientskeep[up,]
patientskeep.down <- patientskeep[down,]

# pull only primary tumour expression

patientskeep.up <- patientskeep.up[,2]
patientskeep.down <- patientskeep.down[,2]

keep <- colnames(tpm.glass) %in% patientskeep.up
keep[1] <- TRUE
keep[2] <- TRUE

tpm.glass.up <- tpm.glass[,keep]

keep <- colnames(tpm.glass) %in% patientskeep.down
keep[1] <- TRUE
keep[2] <- TRUE

tpm.glass.down <- tpm.glass[,keep]

# list of patients

patients.up <- unique(sub(".{3}$","",colnames(tpm.glass.up)[3:ncol(tpm.glass.up)]))
patients.down <- unique(sub(".{3}$","",colnames(tpm.glass.down)[3:ncol(tpm.glass.down)]))

# transcripts of interest

transcripts.int <- tpm.glass[grep("RRBP1",tpm.glass$GeneName),1:2]

# transcripts.int <- read.table("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/batch/FRYL.txt",
#                              header = T)

plots <- vector(mode = "list", length = nrow(transcripts.int))

for (i in 1:nrow(transcripts.int)){
  
  transcript.up <- tpm.glass.up[tpm.glass.up$target_id == transcripts.int[i,1],]
  transcript.up <- as.data.frame(transcript.up)
  head(transcript.up)

  transcript.down <- tpm.glass.down[tpm.glass.down$target_id == transcripts.int[i,1],]
  transcript.down <- as.data.frame(transcript.down)
    
  transcript.df.up <- as.data.frame(matrix(ncol=4, nrow=(ncol(transcript.up)-2)))
  transcript.df.up[,1] <- colnames(transcript.up)[c(-1,-2)]
  colnames(transcript.df.up) <- c("Patient","TPM","Type","Patient.Num")

  for (j in 1:nrow(transcript.df.up)){
    idx <- match(transcript.df.up[j,1],colnames(transcript.up))
    transcript.df.up[j,2] <- transcript.up[1,idx]
    transcript.df.up[j,3] <- "Up"
    transcript.df.up[j,4] <- paste0(j)
  }

  transcript.df.down <- as.data.frame(matrix(ncol=4, nrow=(ncol(transcript.down)-2)))
  transcript.df.down[,1] <- colnames(transcript.down)[c(-1,-2)]
  colnames(transcript.df.down) <- c("Patient","TPM","Type","Patient.Num")

  for (j in 1:nrow(transcript.df.down)){
    idx <- match(transcript.df.down[j,1],colnames(transcript.down))
    transcript.df.down[j,2] <- transcript.down[1,idx]
    transcript.df.down[j,3] <- "Down"
    transcript.df.down[j,4] <- paste0(j+14)
  }


  transcript.df <- rbind(transcript.df.up, transcript.df.down)

  patient <- transcript.df[,4]

# plot bar chart

  p <- ggplot(data = transcript.df, aes(x=Patient.Num, y=TPM)) +
    geom_col(position = "dodge", aes(fill = Type, color = Type)) +
    theme_bw() +
    scale_fill_manual(values = c("Up" = "#190F9733", "Down" = "#208C1733")) +
    scale_color_manual(values = c("Up" = "#190F97", "Down" = "#208C17")) +
    scale_x_discrete(name = "", limits = patient) +
    scale_y_continuous(name = "TPM") +
    ggtitle(transcripts.int[i,2]) +
    theme(axis.text.x = element_text(size = "8"), 
        axis.text.y = element_text(size = "10"),
        axis.title.y = element_text(size = "10", margin = margin(r=5)),
        plot.title = element_text(face = "bold", size=11),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  ggpubr::rremove("grid")

  plots[[i]] <- p
  
}

# SOD2 plots

(plots[[1]] | plots[[2]] | plots[[3]]) /
  (plots[[4]] | plots[[5]] | plots[6]) /
  (plots[[7]] | plots[[8]] | plots[9])

(plots[[10]] | plots[[11]] | plots[[12]]) /
  (plots[[13]])

plots[[19]]


(plots[[1]] | plots[[2]]) 

plots[[1]]

fryl <- plots[[1]]
col1 <- plots[[1]]
col2 <- plots[[2]]
nrxn3 <- plots[[1]]

(col1 | col2) /
  (nrxn3 | fryl )

vtranscripts.int



