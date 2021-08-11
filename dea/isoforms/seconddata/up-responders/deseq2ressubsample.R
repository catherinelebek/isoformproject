library(ggplot2)

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


up <- 6843
down <- 2016


mu <- mean(sig[1:10,1])
sdev <- sd(sig[1:10,1])

tstat <- (mu - 2016) / sqrt((sdev^2)/10)
2*(1 - pt(tstat,9))

1 - pt(tstat,9)

2*pt(-abs(tstat), 9)

quantile <- qt(0.975,9)

CI_lower <- mu - tstat*quantile
CI_upper <- mu + tstat*quantile



p <- ggplot(data = sig, aes(x=rownames(sig),y=Significant)) +
     geom_col(fill = "#190F9733", color = "#190F97") +
     theme_bw() +
     scale_x_discrete(name = "Iteration", limits = c("1","2","3","4","5","6","7","8","9","10")) +
     scale_y_continuous(name = "Significant Isoforms", limits = c(0,7000)) +
     theme(axis.text.x = element_text(size = "12"), 
           axis.text.y = element_text(size = "12"),
           axis.title.x = element_text(face = "bold", size = "12", margin = margin(t=20)),
           axis.title.y = element_text(face = "bold", size = "12", margin = margin(r=20)),
           plot.margin = unit(c(1,1,0.5,0.5), "cm")) +
     geom_hline(yintercept = up, color = "black", size = 1, linetype = "dotdash") +
     geom_hline(yintercept = down, color = "black", size = 1, linetype = "dashed") +
     geom_hline(yintercept = mu, color = "red", size = 1.2, line) +
     ggpubr::rremove("grid")

p

mu




   