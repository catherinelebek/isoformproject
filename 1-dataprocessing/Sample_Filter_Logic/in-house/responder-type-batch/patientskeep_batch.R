# import metadata

metadata <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/MetaData_LS_230621.txt",
                     header = TRUE, "\t")

# import sequencing data

sequencing <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/SequencingMetrics_LS_230621.txt",
                       header = TRUE, "\t")

# add patient column to sequencing data


sequencing$Patient.ID <- sub("rimary","",sequencing$Sample)
sequencing$Patient.ID <- sub("ecurrent","",sequencing$Patient.ID)
sequencing$Patient.ID <- gsub(".{2}$","",sequencing$Patient.ID)

# remove patients that are not in both sequencing and metadata (expect this to be 6 samples)

metadata.mod <- metadata[metadata$Patient.ID %in% sequencing$Patient.ID,]
sequencing.mod <- sequencing[sequencing$Patient.ID %in% metadata$Patient.ID,]

# Logic for patients that should be removed
# Initial diagnosis of primary GBM
# Same location of primary and recurrent tumour, and not NA or empty
# Radiotherapy and TMZ
# IDH 0

keep <- metadata.mod$Initial.Diagnosis == "Primary GBM" & metadata.mod$Location.Primary == metadata.mod$Location.Recurrence &
  !is.na(metadata.mod$Location.Primary) & metadata.mod$Location.Primary != "TBC" & 
  (metadata.mod$Non.Surgical.Treatment == "Radiotherapy and TMZ" | metadata.mod$Non.Surgical.Treatment == "Radiotherapy and TMZ +") &
  metadata.mod$IDH == 0

patientskeep.metadata <- metadata.mod[keep,3]
patientskeep.metadata.P <- paste(patientskeep.metadata,"_P", sep = "")
patientskeep.metadata.R <- paste(patientskeep.metadata,"_R", sep = "")
patientskeep.metadata <- c(patientskeep.metadata.P, patientskeep.metadata.R)

# get list of patients with read below above 30m

keep <- sequencing.mod$Reads < 30000000
sequencing.mod <- sequencing.mod[keep,1]
sequencing.mod <- sub("rimary","",sequencing.mod)
sequencing.mod <- sub("ecurrent", "", sequencing.mod)

# remove any patients from patientskeep.metadata that have at least one of paired samples below 30m reads

keep <- !patientskeep.metadata %in% sequencing.mod
patientskeep <- patientskeep.metadata[keep]

# finally, check that sequencing data exists for both samples for 66 patients in patientskeep

temp <- c()

for (i in 1:length(patientskeep)){
  temp[i] <- sum(sequencing$Sample == patientskeep[i])
}

table(temp)

# write list of patients to keep to file

write.table(patientskeep, 
            "/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/patientskeep_batch.txt",
            col.names = F, row.names = F)
