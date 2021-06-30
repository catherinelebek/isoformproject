library(DESeq2)

# load DESeq2 object as dds

load("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2.RData")

# load full list of transcripts in order to pull through gene names

genelist <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

# extract just transcripts and gene names

genelist <- genelist[,c(1,10)]

# create results table from dds object

res <- results(dds, alpha = 0.05)

# order results by raw p-value

resOrdered <- res[order(res$pvalue),]

# put in a dataframe so can merge with the list of transcripts and gene names

resOrdered <- as.data.frame(resOrdered)

# run merge to pull through gene names into the DESeq2 results table

merge <- merge(resOrdered, genelist, by.x = "row.names", by.y = "EnsID", all.x = TRUE)

# again order by p-value

merge <- merge[order(res$pvalue),]

# move gene list to be second column

merge <- merge[,c(1,8,2:7)]

# get a summary of the results at the 90% confidence level

summary(res)

# save results to csv

write.csv(merge, "~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/deseq2results.csv")

