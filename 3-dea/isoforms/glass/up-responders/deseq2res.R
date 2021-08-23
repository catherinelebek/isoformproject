# This script processes the results of the paired-sample DEA for up-responders (glass dataset)

library(DESeq2)

# load DESeq2 object as dds

load("~/Documents/Semester3/Project/Results/dea/isoforms/glass/up-responders/glassfilter/deseq2.RData")

# create results table from dds object

res <- results(dds, alpha = 0.05)

# order results by raw p-value

resOrdered <- res[order(res$pvalue),]

# put in a dataframe so can merge with the list of transcripts and gene names

resOrdered <- as.data.frame(resOrdered)

# get a summary of the results at the 90% confidence level

summary(res)

# save results to csv

write.csv(resOrdered, "/Users/catherinehogg/Documents/Semester3/Project/Results/dea/isoforms/glass/up-responders/glassfilter/deseq2results.csv")


table(resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) > 1)


# Are any isoforms of YOD1 significant?
resOrdered[rownames(resOrdered) == "ENST00000315927",] # YOD1-201
resOrdered[rownames(resOrdered) == "ENST00000367084",] # YOD1-202

