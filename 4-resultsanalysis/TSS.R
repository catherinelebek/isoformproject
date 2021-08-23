library(GSA)
library(dplyr)

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



head(transcripts.tss) # transcript IDs annotated with TSS
head(jarid.tss) # jarid2 TSS
head(jarid.tss.transcripts) # transcript IDs with JARID2 TSS
jarid2.gene <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                          header = F, sep = ",") # JARID2 gene IDs

jarid2.gene <- c(t(jarid2.gene))

mapping <- read.delim("~/Documents/Semester3/Project/InputData/output.txt", header = F, sep = "\t") # gene-to-transcript IDs

colnames(mapping) <- c("Gene.EnsID","Transcript.EnsID")

head(mapping)

mapping$transcript.inc <- mapping$Transcript.EnsID %in% transcripts.tss$Transcript
mapping$transcript.jarid.tss <- mapping$Transcript.EnsID %in% jarid.tss.transcripts
table(mapping$transcript.jarid.tss)
mapping$transcript.gene <- sub("\\..*","",mapping$Gene.EnsID) %in% jarid2.gene
table(mapping$transcript.gene)

test <- mapping[mapping$transcript.gene == TRUE,1]
test <- sub("\\..*","",test)
length(unique(test))

mapping.count <- mapping %>% group_by(Gene.EnsID) %>% summarise(transcript.inc = n(), jarid2gene = sum(transcript.gene), 
                                                                                      jarid2tss = sum(transcript.jarid.tss))

mapping.count

nrow(mapping.count[mapping.count$jarid2gene > mapping.count$jarid2tss,])
nrow(mapping.count[mapping.count$jarid2gene > 0,])
nrow(mapping.count[mapping.count$jarid2gene == mapping.count$jarid2tss & mapping.count$jarid2gene != 0,])
nrow(mapping.count[mapping.count$jarid2gene > 0 & mapping.count$jarid2tss == 0,])

2953/5234

length(unique(mapping.count$Gene.EnsID))
