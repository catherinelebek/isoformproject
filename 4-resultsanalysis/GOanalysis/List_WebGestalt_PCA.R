# input

input <- "topisoforms_glass.csv"

# all isoforms

output.all <- "pc4.all.txt"

# top isoforms

output.top.up <- "pca4.top.up.txt"
output.top.down <- "pca4.top.down.txt"

# make list

list.interest <- read.delim(paste0("~/Documents/Semester3/Project/Results/resultsanalysis/PCA/",input),
                            sep = ",", header = T)

# keep only the first 4 PCs

list.interest <- list.interest[,c(1,2,3,4,5)]

# load full list of transcripts in order to pull through gene names

genelist <- read.delim("/Users/catherinehogg/Documents/Semester3/Project/InputData/isoforms/seconddata/PvR_isoformCounts_all_LS_23062021.txt.txt",
                       header = T, sep = "\t")

# extract just transcripts and gene names

genelist <- genelist[,c(1,10)]

genelist$EnsIDSimp <- sub("\\..*","",genelist$EnsID)

table(list.interest$X %in% genelist$EnsIDSimp)

# run a merge

list.interest <- merge(list.interest, genelist, by.x = "X", by.y = "EnsIDSimp", all.x = TRUE)

# get simple gene name

list.interest$GeneName <- sub("-.*","",list.interest$GeneName)

# order based on actual values

list.interest <- list.interest[order(list.interest$PC3, decreasing = T),]

head(list.interest)

cols <- list.interest[is.na(list.interest$GeneName),1]

# find most important isoforms

plot(list.interest$PC3)
qqnorm(list.interest$PC3)
qqline(list.interest$PC3)
abline(a = 0.008, b = 0, col = "red")
abline(a = -0.0085, b = 0, col = "red")


#PC3 = -0.0085, 0.006
#PC4 = -0.0085, 0.006


# get list of most important isoforms

list.interest.top.up <- list.interest[list.interest$PC3 > 0.008,7]
head(list.interest)

list.interest.top.down <- list.interest[list.interest$PC3 < -0.0085,7]


# create list of unique values

list.interest.all <- unique(list.interest$GeneName)
list.interest.top.up <- unique(list.interest.top.up)
list.interest.top.down <- unique(list.interest.top.down)

# remove NA

list.interest.all <- list.interest.all[!is.na(list.interest.all)]
list.interest.top.up <- list.interest.top.up[!is.na(list.interest.top.up)]
list.interest.top.down <- list.interest.top.down[!is.na(list.interest.top.down)]


# write to file

write.table(list.interest.all, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_PCA/",output.all),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.top.up, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_PCA/",output.top.up),
            sep = "\t", row.names = F, col.names = F, quote = F)

write.table(list.interest.top.down, paste0("~/Documents/Semester3/Project/Results/resultsanalysis/GO_PCA/",output.top.down),
            sep = "\t", row.names = F, col.names = F, quote = F)
