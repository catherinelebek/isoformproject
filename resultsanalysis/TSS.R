library(GSA)

jarid.tss <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/TSS.gmt",
                        sep = "\t", header = T)

transcripts.tss <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/transcript_to_tss_position.txt",
                        sep = "\t", header = F)

jarid.tss <- jarid.tss[-1,]
colnames(transcripts.tss) <- c("Transcript","TSS")
  
head(transcripts.tss)

jarid.tss.transcripts <- transcripts.tss[transcripts.tss$TSS %in% jarid.tss,1]

write.csv(jarid.tss.transcripts, "/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/jarid.tss.transcripts.csv",
          row.names = F)
