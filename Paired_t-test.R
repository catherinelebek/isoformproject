# import count data

dat <- read.delim("PvR_isoformCounts_all.txt", header = TRUE)
dat_full <- dat
head(dat_full)

# rearrange columns

dat <- dat_full[,c(-10,-11)]
rownames(dat) <- dat[,1]
dat <- dat[,c(-1)]
head(dat)

# remove unwanted samples

patients.remove <- read.delim("patients_remove.txt", header = F)

# convert list to vector

patients.remove <- as.vector(t(patients.remove))

# filter

keep <- !sub("_.*","",colnames(dat)) %in% patients.remove
dat <- dat[,keep]

# remove rows that add up to zero

dat <- dat[rowSums(dat[,1:ncol(dat)]) != 0,]

# split data into primary and recurrent

primary <- sub("_R.*","",colnames(dat)) %in% colnames(dat)
primary <- dat[,primary]
colnames(primary) <- sub("_.*","",colnames(primary))

recurrent <- sub("*_P.*","",colnames(dat)) %in% colnames(dat)
recurrent <- dat[,recurrent]
colnames(recurrent) <- sub("_.*","",colnames(recurrent))

# column names should be equal

table(colnames(primary) == colnames(recurrent))

# create a list with one item for every transcript
# each transcript has a data frame of counts for each patient's primary and recurrent tumour
# expression delta calculated for each patient and transcript
# this will enable me to run a paired t-test for each transcript

dat_inc <- nrow(dat)

dat_list <- vector(mode = "list", length = dat_inc) # create empty list
names(dat_list) <- rownames(dat[1:dat_inc,]) # name list elements as transcripts

# add empty data frame to each list element, with a column for primary, recurrent, delta and LFC expression
# and a row for every patient

for (i in 1:length(dat_list)){
  temp <- data.frame(matrix(ncol=4,nrow=ncol(primary)))
  rownames(temp) <- colnames(primary[1:ncol(primary)])
  colnames(temp) <- c("P","R","D","LFC")
  dat_list[[i]] <- temp
}

# add data to the data frames in the list

for (i in 1:length(dat_list)){
  for (j in 1:ncol(primary)){
  dat_list[[i]][j,1] <- primary[i,j] # primary tumour expression
  dat_list[[i]][j,2] <- recurrent[i,j] # recurrent tumour expression
  dat_list[[i]][j,3] <- (dat_list[[i]][j,2] - dat_list[[i]][j,1]) # delta expression
  dat_list[[i]][j,4] <- log((dat_list[[i]][j,2])/(dat_list[[i]][j,1]),2) # LFC expression
 }
}

# calculate correlations for each transcript

corr <- rep(NA,length(dat_list)) 

for (i in 1:length(dat_list)){
  corr[i] <- cor(dat_list[[i]][,1],dat_list[[i]][,2])
}

write(corr, "correlations_per_patient.txt")

# do a paired t-test for each transcript
# output values into dataframe

ttest <- data.frame(matrix(ncol=8,nrow=length(dat_list)))
colnames(ttest) <- c("p.value","LFC","dof","statistic","CI-lower","CI-upper","Estimate","Sterr")
rownames(ttest) <- names(dat_list)

# output results from ttest into a data frame

for (i in 1:length(dat_list)){
  temp <- t.test(dat_list[[i]][,3], alternative = "two.sided", mu = 0, conf.level = 0.95)
  ttest[i,1] <- temp$p.value
  ttest[i,2] <- mean(dat_list[[i]][,4])
  ttest[i,3] <- temp$parameter
  ttest[i,4] <- temp$statistic
  ttest[i,5] <- temp$conf.int[1]
  ttest[i,6] <- temp$conf.int[2]
  ttest[i,7] <- temp$estimate
  ttest[i,8] <- temp$stderr
}

adj.pval <- p.adjust(ttest$p.value,"BH")
ttest$adj.pval <- adj.pval # multiplicity adjustment

# merge with data frame containing gene name and gene type

ttest <- merge(ttest, dat_full[,c("EnsID","GeneName","GeneType")], by.x = "row.names", by.y = "EnsID")
rownames(ttest) <- ttest[,1]
ttest <- ttest[,-1]
ttest <- ttest[,c(10,11,1,9,2:8)]

# reorder by increasing p-value

ttest <- ttest[order(ttest$adj.pval),]
write.table(ttest, "paired_ttest_results.csv")
