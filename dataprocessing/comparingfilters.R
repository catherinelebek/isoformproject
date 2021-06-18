# read in data from georgette's scripts

gtinc <- gt
gt <- datfull[!datfull$EnsID %in% gtinc, 1]

ch.counts <- read.csv("/Users/catherinehogg/Documents/Semester3/Project/Results/localresults/lowexpressionomit.csv", header = TRUE, sep = " ")
ch.counts <- as.vector(t(ch.counts))

table(gt %in% ch.counts)
table(ch.counts %in% gt)

idx <- (!ch.counts %in% gt)
err <- ch.counts[idx]

test <- err[1]
test


ynorm[y$genes$EnsID == test,] 
prim <- ynormprimary[y$genes$EnsID == test,]
recc <- ynormrecurrent[y$genes$EnsID == test,]

primcount <- sum(prim > lowerq)
recccout <- sum(recc > lowerq)


head(dat.fpkm)

prim1 <- datfpkmprimary[rownames(datfpkmprimary) == test,]
recc1 <- datfpkmrecurrent[rownames(datfpkmrecurrent) == test,]

primcount1 <- sum(prim1 > lowerq)
recccout1 <- sum(recc1 > lowerq)


