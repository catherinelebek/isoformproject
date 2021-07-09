results <- read.delim("~/Documents/Semester3/Project/Results/resultsanalysis/summarydf.csv", header = T, sep = ",")
head(results)

jarid2 <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                     header = F, sep = ",")
jarid2 <- c(t(jarid2))

isoforms <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv", header = T, sep = ",")
head(isoforms)

isoforms.up <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv", header = T, sep = ",")
  
isoforms.down <- read.delim("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/down-responders/deseq2results.csv", header = T, sep = ",")

# potentially interesting genes/transcripts

head(results)

#SCN8A

results[results$GeneName == "DCX",]
isoforms.all[sub("-.*","",isoforms.all$GeneName) == "OPALIN",]
isoforms.up[sub("-.*","",isoforms.up$GeneName) == "PTPRD",]
isoforms.down[sub("-.*","",isoforms.down$GeneName) == "PTPRD",]

wang <- c("MAP4K4","MAPK8","SRSF5","PDGFRB","PIK3R1","FOS","MAPK9","TGFB1","SMAD2","SMAD4",
          "NFKB1","RELA","STAT1","STAT3","PTPRZ1","CD44","MAP4K4","MAPK9","MAPK10",
          "FGFR1","FGFR2","EGFR","TNC","FN1","SRSF5","SRSF9")

wang <- c("MAP4K4","MAPK9","MAPK10",
          "FGFR1","FGFR2","EGFR","TNC","FN1","SRSF5","SRSF9")
results[results$GeneName %in% wang,]

head(results)

test <- read.delim("PvR_isoformCounts_all_LS_23062021.txt.txt", header = T,
                   sep = "\t")
head(test)

