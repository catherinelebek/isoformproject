library(DESeq2)
# library(BiocParallel)

# counts.full <- read.delim("~/Documents/Semester3/Project/Results/filtered_data/isoforms/seconddata/up-responders/PvR_isoformCounts_filtered.txt",header = T, sep = "\t")
counts.full <- read.delim("/nobackup/bs20chlb/inputdata/seconddata/up-responders/PvR_isoformCounts_filtered.txt",header = T, sep = "\t")

counts <- counts.full

rownames(counts) <- counts[,1]
counts <- counts[,c(-1,-2,-3)]

head(counts)

samples <- data.frame(matrix(ncol = 2, nrow = ncol(counts)))
colnames(samples) <- c("patientid","tumourtype")
rownames(samples) <- colnames(counts)

for (i in 1:nrow(samples)){
  samples[i,1] <- gsub(".{2}$","",colnames(counts)[i])
  samples[i,2] <- substr(colnames(counts)[i], nchar(colnames(counts)[i]), nchar(colnames(counts)[i]))
}

samples

all(rownames(samples) == colnames(counts))

# define how many iterations to do

n <- 10

# loop deseq2 analysis over sub-sample iterations

for (i in 1:n){

  cols <- unique(sub(".{2}$","",rownames(samples))) 
  idx <- sample(1:length(cols),22)
  cols <- cols[idx]

  counts.sub <- counts[,samples$patientid %in% cols]
  samples.sub <- samples[samples$patientid %in% cols,]

  samples.sub$patientid <- as.factor(samples.sub$patientid)
  samples.sub$tumourtype <- as.factor(samples.sub$tumourtype)

  all(rownames(samples.sub) == colnames(counts.sub))

  # create DESeq data set

  dds <- DESeqDataSetFromMatrix(countData = counts.sub,
                              colData = samples.sub,
                              design = ~ patientid + tumourtype)


  dds$tumourtype <- relevel(dds$tumourtype, ref = "P")

 # dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(4))

  dds <- DESeq(dds)
  
  # create results table from dds object
  
  res <- results(dds, alpha = 0.05)
  
  # order results by raw p-value
  
  resOrdered <- res[order(res$pvalue),]
  
  # put in a dataframe so can merge with the list of transcripts and gene names
  
  resOrdered <- as.data.frame(resOrdered)
  
  # save results to csv

  write.csv(resOrdered, paste0("/nobackup/bs20chlb/outputdata/dea/seconddata/up-responders/subsample/deseq2results",i,".csv"))
  
}
  
  
  
