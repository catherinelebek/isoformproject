check <- read.delim("~/Documents/Semester3/Project/Results/localresults/filter3/PvR_isoformCounts_filtered.txt",header = T, sep = " ")

transcripts <- vector(mode = "list", length = nrow(check))

# create an empty dataframe for every transcript

for (i in 1:length(transcripts)){
  temp <- data.frame(matrix(ncol = 2, nrow = 43))
  colnames(temp) <- c("Primary","Recurrent")
  rownames(temp) <- unique(sub("_.*","",colnames(check[,4:ncol(check)])))
  transcripts[[i]] <- temp
}


for (i in 1:length(transcripts)){
  for (j in 1:nrow(transcripts[[1]])){
  transcripts[[i]][j,1] <- check[check$EnsID == names(transcripts)[i], colnames(check) == paste(rownames(transcripts[[i]][j,]),"_P",sep="")]
  transcripts[[i]][j,2] <- check[check$EnsID == names(transcripts)[i], colnames(check) == paste(rownames(transcripts[[i]][j,]),"_R",sep="")]
  }
}