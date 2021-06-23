fullcounts <- read.delim("/nobackup/bs20chlb/inputdata/archive/PvR_isoformCounts_all.txt",header = T, sep = " \t")

transcripts <- vector(mode = "list", length = nrow(fullcounts))

# create an empty dataframe for every transcript

for (i in 1:length(transcripts)){
  temp <- data.frame(matrix(ncol = 2, nrow = 43))
  colnames(temp) <- c("Primary","Recurrent")
  rownames(temp) <- unique(sub("_.*","",colnames(fullcounts[,4:ncol(fullcounts)])))
  transcripts[[i]] <- temp
}


for (i in 1:length(transcripts)){
  for (j in 1:nrow(transcripts[[1]])){
  transcripts[[i]][j,1] <- fullcounts[fullcounts$EnsID == names(transcripts)[i], colnames(fullcounts) == paste(rownames(transcripts[[i]][j,]),"_P",sep="")]
  transcripts[[i]][j,2] <- fullcounts[fullcounts$EnsID == names(transcripts)[i], colnames(fullcounts) == paste(rownames(transcripts[[i]][j,]),"_R",sep="")]
  }
}
