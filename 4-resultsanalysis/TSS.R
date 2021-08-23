# This script has two purposes:
# 1. Produce a list of transcripts with JARID2 TSSs
# 2. Calculate % of JBSgenes that have both transcripts with JARID2 TSSs and without

library(dplyr)

# read in list of JARID2 TSSs
jarid.tss <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/TSS.gmt",
                        sep = "\t", header = T)

jarid.tss <- jarid.tss[-1,] # remove NA at top of list

# read in list of transcripts and corresponding TSSs
transcripts.tss <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/transcript_to_tss_position.txt",
                        sep = "\t", header = F)

colnames(transcripts.tss) <- c("Transcript","TSS")
  
head(transcripts.tss)

# find transcripts with JARID2 TSSs and write to csv

jarid.tss.transcripts <- transcripts.tss[transcripts.tss$TSS %in% jarid.tss,1]

write.csv(jarid.tss.transcripts, "/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/jarid.tss.transcripts.csv",
          row.names = F)



# read in list of JBSgenes

jarid2.gene <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/genes/jarid2.csv",
                          header = F, sep = ",") # JARID2 gene IDs

jarid2.gene <- c(t(jarid2.gene))

# read in mapping between gene and transcript IDs

mapping <- read.delim("~/Documents/Semester3/Project/InputData/output.txt", header = F, sep = "\t") # gene-to-transcript IDs

colnames(mapping) <- c("Gene.EnsID","Transcript.EnsID")

head(mapping)

# add column to mapping dataframe to indicate if transcript has a JARID2 TSS

mapping$transcript.jarid.tss <- mapping$Transcript.EnsID %in% jarid.tss.transcripts
table(mapping$transcript.jarid.tss)

# add columnn to mapping dataframe to indicate if gene relating to transcript is a JBSgene

mapping$transcript.gene <- sub("\\..*","",mapping$Gene.EnsID) %in% jarid2.gene
table(mapping$transcript.gene)

# filter mapping dataframe to only include transcripts that are from JBSgenes, and extract the gene IDs

test <- mapping[mapping$transcript.gene == TRUE,1]
test <- sub("\\..*","",test)

# extract unique list of JBSgenes to check correct number

length(unique(test))

# use dplyr package to count how many are genes are JBSgenes, how many transcripts relate to each gene, and how many of these have JARID TSSs

mapping.count <- mapping %>% group_by(Gene.EnsID) %>% summarise(transcript.inc = n(), jarid2gene = sum(transcript.gene), 
                                                                                      jarid2tss = sum(transcript.jarid.tss))

mapping.count

nrow(mapping.count[mapping.count$jarid2gene > mapping.count$jarid2tss,]) # count of JBSgenes that have both transcripts with JARID2 TSSs and without
nrow(mapping.count[mapping.count$jarid2gene > 0,]) # count of JBSgenes
nrow(mapping.count[mapping.count$jarid2gene == mapping.count$jarid2tss & mapping.count$jarid2gene != 0,]) # count of JBSgenes only with transcripts with JARID2 TSSs
nrow(mapping.count[mapping.count$jarid2gene > 0 & mapping.count$jarid2tss == 0,]) # counts of JBSgenes with no transcripts with JARID2 TSSs

2953/5234 # percentage of JBSgenes with two types of transcripts: those with JARID2 TSSs and those without

