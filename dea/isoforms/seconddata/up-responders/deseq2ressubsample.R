sig <- as.data.frame(matrix(ncol = 3, nrow = 10))
colnames(sig) <- c("Significant", "Up", "Down")

for (i in 1:10){

 results <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/subsample/deseq2results",i,".csv"), 
                       header = T, sep = ",")
 sig[i,1] <- sum(results$padj < 0.05, na.rm = T)
 sig[i,2] <- sum(results$padj < 0.05 & results$log2FoldChange > 0, na.rm = T)
 sig[i,3] <- sum(results$padj < 0.05 & results$log2FoldChange < 0, na.rm = T)

}

sig$percentup <- sig$Up/sig$Significant


sig[11,] <- c(6843,4692,2151,4692/6842)
sig[12,] <- c(2020,246,1774, 246/2020)

sig$type <- c(rep("Subsample",10),"All up-responders", "All down-responders")
sig$type <- as.factor(sig$type)

sig

rownames(sig) <- c("1","2","3","4","5",
                   "6","7","8","9","10",
                   "U","D")

barplot(sig$Significant, ylim = c(0,7000), col = "lightblue", pch = 16, names.arg = rownames(sig))

mu <- mean(sig[1:10,1])
sdev <- sd(sig[1:10,1])

tstat <- (mu - 2020) / sqrt((sdev^2)/10)
2*(1 - pt(tstat,9))

sig$percentup <- sig$Up/sig$Significant
sig


2*pt(-abs(tstat), 9)

quantile <- qt(0.975,9)

CI_lower <- mu - tstat*quantile
CI_upper <- mu + tstat*quantile


barplot(sig$Significant, ylim = c(0,7000), col = "lightblue", pch = 16, names.arg = rownames(sig))
abline(mu,0, col = "red")
abline(CI_lower, 0, col = "red")
abline(CI_upper, 0, col = "red")
abline(2020,0,col = "green")
   