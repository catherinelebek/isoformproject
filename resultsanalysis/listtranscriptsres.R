load("~/Documents/Semester3/Project/Results/resultsanalysis/transcripts.RData")

# This list has every transcript for filter1 with the primary and recurrent counts
# There are 43 patients for each transcripts, both with primary and recurrent samples

transcripts$ENST00000085219.9

missing

dat <- datfull
test <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/Scripts/isoformproject/local/localdata/PvR_isoformCounts_all.txt",
           header = T,
           sep = "\t")

test[test$EnsID == "ENST00000552351.5",]

dat <- dat["HIF1A" == sub("-.*","",dat$GeneName),]
dat


