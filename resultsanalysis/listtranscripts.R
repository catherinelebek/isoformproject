check <- read.delim("/nobackup/bs20chlb/inputdata/filter1/PvR_isoformCounts_filtered.txt",header = T, sep = " ")

transcripts <- vector(mode = "list", length = nrow(check))
names(transcripts) <- check[,1]

func <- function(x){
    temp <- which(names(transcripts) == x)
}


# create an empty dataframe for every transcript

for (i in 1:length(transcripts)){
  temp <- data.frame(matrix(ncol = 2, nrow = (ncol(check[,4:ncol(check)]))/2))
  colnames(temp) <- c("Primary","Recurrent")
  rownames(temp) <- unique(sub("_.*","",colnames(check[,4:ncol(check)])))
  transcripts[[i]] <- temp
}

num <- lapply(names(transcripts),func)

temp2 <- function(row, col, x){
  row <- as.numeric(rownames(x)
  col <- as.numeric(colnames(x)
  temp <- check$EnsID == names(transcripts)[num], colnames(check) == paste(rownames(transcripts)[[num]]
}
                                                                 
for (i in 1:length(transcripts)){
  for (j in 1:nrow(transcripts[[i]])){
  transcripts[[i]][j,1] <- check[check$EnsID == names(transcripts)[i], colnames(check) == paste(rownames(transcripts[[i]][j,]),"_P",sep="")]
  transcripts[[i]][j,2] <- check[check$EnsID == names(transcripts)[i], colnames(check) == paste(rownames(transcripts[[i]][j,]),"_R",sep="")]
  }
}

save(transcripts, file = "transcripts.RData")

