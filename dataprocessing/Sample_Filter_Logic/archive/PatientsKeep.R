# import metadata

metadata <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/firstdata/samplefilters/MetaData.csv",
                     header = TRUE, ",")

# import sequencing data

sequencing <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/firstdata/samplefilters/SequencingMetrics.txt",
                     header = TRUE, "\t")


# add patient column to sequencing data

sequencing$Patient.ID <- sub("rimary","",sequencing$Sample)
sequencing$Patient.ID <- sub("ecurrent","",sequencing$Patient.ID)
sequencing$Patient.ID <- gsub(".{2}$","",sequencing$Patient.ID)


# remove patients that are not in both sequencing and metadata

metadata.mod <- metadata[metadata$Patient.ID %in% sequencing$Patient.ID,]
sequencing.mod <- sequencing[sequencing$Patient.ID %in% metadata$Patient.ID,]

# Logic for patients that should be removed
# Inital diagnosis of primary GBM
# Same location of primary and recurrent tumour, and not NA or empty
# Radiotherapy and TMZ
# IDH 0

keep <- metadata.mod$Initial.Diagnosis == "Primary GBM" & metadata.mod$Location.Primary == metadata.mod$Location.Recurrence &
        !is.na(metadata.mod$Location.Primary) & metadata.mod$Location.Primary != "TBC" & 
        (metadata.mod$Non.Surgical.Treatment == "Radiotherapy and TMZ" | metadata.mod$Non.Surgical.Treatment == "Radiotherapy and TMZ +") &
        metadata.mod$IDH == 0

patientskeep.metadata <- metadata.mod[keep,3]

# get list of patients with read below above 30m

keep <- sequencing.mod$Reads < 30000000
sequencing.mod <- sequencing.mod[keep,17]
patientsremove.below30 <- unique(sequencing.mod) # because if one paired sample below 30m reads then need to remove both

# keep any patients from patientskeep.metadata that have at both paired samples above 30m reads

keep <- !patientskeep.metadata %in% patientsremove.below30
patientskeep <- patientskeep.metadata[keep]

# finally, check that sequencing data exists for both samples for 66 patients in patientskeep

temp <- c()

for (i in 1:length(patientskeep)){
  temp[i] <- sum(sequencing$Patient.ID == patientskeep[i])
}

table(temp)

# write list of patients to keep to file

write.table(patientskeep, 
            "/Users/catherinehogg/Documents/Semester3/Project/Results/filtered_data/isoforms/firstdata/patientskeep.txt",
            col.names = F, row.names = F)
