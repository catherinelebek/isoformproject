# This script is to select and output genes for use in BP GO term enrichment analysis based on results from paired-sample isoform-level DEAs

# first for up-responders #####

input <- "up-responders/deseq2results_merge.csv"

# all isoforms

output.all <- "up.isoforms.all.txt"

# all isoforms above threshold

output.sig.up <- "up.isoforms.sig.up.txt"
output.sig.down <- "up.isoforms.sig.down.txt"

# all isoforms above threshold and JARID2 associated, but without a jarid2 tss

output.sig.up.jarid.ass <- "up.isoforms.sig.up.jarid.ass.txt"
output.sig.down.jarid.ass <- "up.isoforms.sig.down.jarid.ass.txt"

# all isoforms above threshold and JARID2 TSS associated

output.sig.up.jarid.tss <- "up.isoforms.sig.up.jarid.tss.txt"
output.sig.down.jarid.tss  <- "up.isoforms.sig.down.jarid.tss.txt"

# all isoforms above threshold and not JARID2 associated

output.sig.up.notjarid <- "up.isoforms.sig.up.notjarid.txt"
output.sig.down.notjarid  <- "up.isoforms.sig.down.notjarid.txt"


# now for down-responders ####

input <- "down-responders/deseq2results_merge.csv"

# all isoforms

output.all <- "down.isoforms.all.txt"

# all isoforms above threshold

output.sig.up <- "down.isoforms.sig.up.txt"
output.sig.down <- "down.isoforms.sig.down.txt"

# all isoforms above threshold and JARID2 associated

output.sig.up.jarid.ass <- "down.isoforms.sig.up.jarid.ass.txt"
output.sig.down.jarid.ass <- "down.isoforms.sig.down.jarid.ass.txt"

# all isoforms above threshold and JARID2 TSS associated

output.sig.up.jarid.tss <- "down.isoforms.sig.up.jarid.tss.txt"
output.sig.down.jarid.tss  <- "down.isoforms.sig.down.jarid.tss.txt"

# all isoforms above threshold and not JARID2 associated

output.sig.up.notjarid <- "down.isoforms.sig.up.notjarid.txt"
output.sig.down.notjarid  <- "down.isoforms.sig.down.notjarid.txt"


#####

list.interest <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/",input),
                          sep = ",", header = T)

list.interest$GeneSymbol <- sub("-.*","",list.interest$GeneName)
head(list.interest)

list.interest.sig.up <- list.interest[list.interest$log2FoldChange > 1 &
                                   list.interest$padj < 0.05,18]

list.interest.sig.down <- list.interest[list.interest$log2FoldChange < -1 &
                                        list.interest$padj < 0.05,18]

list.interest.sig.up.jarid.ass <- list.interest[list.interest$log2FoldChange > 1 &
                                                  list.interest$padj < 0.05 &
                                                  (list.interest$jarid2.gene == TRUE & list.interest$jarid2.tss == FALSE),18]

list.interest.sig.down.jarid.ass <- list.interest[list.interest$log2FoldChange < -1 &
                                                    list.interest$padj < 0.05 &
                                                    (list.interest$jarid2.gene == TRUE & list.interest$jarid2.tss == FALSE),18]


list.interest.sig.up.jarid.tss <- list.interest[list.interest$log2FoldChange > 1 &
                                                 list.interest$padj < 0.05 &
                                                 list.interest$jarid2.tss == TRUE,18]


list.interest.sig.down.jarid.tss <- list.interest[list.interest$log2FoldChange < -1 &
                                                    list.interest$padj < 0.05 &
                                                    list.interest$jarid2.tss == TRUE,18]

list.interest.sig.up.notjarid <- list.interest[list.interest$log2FoldChange > 1 &
                                                  list.interest$padj < 0.05 &
                                                  list.interest$jarid2.tss == FALSE &
                                                  list.interest$jarid2.gene == FALSE,18]


list.interest.sig.down.notjarid <- list.interest[list.interest$log2FoldChange < -1 &
                                                 list.interest$padj < 0.05 &
                                                 list.interest$jarid2.tss == FALSE &
                                                 list.interest$jarid2.gene == FALSE,18]



list.interest.sig.up <- unique(list.interest.sig.up)
list.interest.sig.down <- unique(list.interest.sig.down)
list.interest.sig.up.jarid.ass <- unique(list.interest.sig.up.jarid.ass)
list.interest.sig.down.jarid.ass <- unique(list.interest.sig.down.jarid.ass)
list.interest.sig.up.jarid.tss <- unique(list.interest.sig.up.jarid.tss)
list.interest.sig.down.jarid.tss <- unique(list.interest.sig.down.jarid.tss)
list.interest.sig.up.notjarid <- unique(list.interest.sig.up.notjarid)
list.interest.sig.down.notjarid <- unique(list.interest.sig.down.notjarid)


list.interest.all <- list.interest$GeneSymbol
list.interest.all <- unique(list.interest.all)


write.table(list.interest.sig.up, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.sig.up),
          sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.sig.down, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.sig.down),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.sig.up.jarid.ass, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.sig.up.jarid.ass),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.sig.down.jarid.ass, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.sig.down.jarid.ass),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.sig.up.jarid.tss, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.sig.up.jarid.tss),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.sig.down.jarid.tss, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.sig.down.jarid.tss),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.sig.up.notjarid, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.sig.up.notjarid),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.sig.down.notjarid, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.sig.down.notjarid),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.all, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/",output.all),
            sep = "\t", row.names = F, col.names = F, quote = F)


# for all patients


list.interest <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/deseq2results.csv"),
                            sep = ",", header = T)
list.interest$GeneIDSimp <- sub("-.*","",list.interest$GeneName)
list.interest.all <- list.interest$GeneIDSimp
list.interest.all <- unique(list.interest.all)
write.table(list.interest.all, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/all.isoforms.all.txt"),
            sep = "\t", row.names = F, col.names = F, quote = F)

list.interest.justall <-  read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/justall.txt"),
                                     sep = "\t", header = F)

list.interest.justall <- merge(list.interest.justall, list.interest[,c("Row.names","GeneIDSimp")],
                               by.x = "V1", by.y = "Row.names", all.x = TRUE)
list.interest.justall <- list.interest.justall[,2]

write.table(list.interest.justall, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/all.isoforms.justall.txt"),
            sep = "\t", row.names = F, col.names = F, quote = F)
