library(ggplot2)

# create empty data frame to store results from each iteration

sig <- as.data.frame(matrix(ncol = 3, nrow = 10))
colnames(sig) <- c("Significant", "Up", "Down")

# import results

for (i in 1:10){

 results <- read.delim(paste0("~/Documents/Semester3/Project/Results/dea/isoforms/seconddata/up-responders/subsample/deseq2results",i,".csv"), 
                       header = T, sep = ",")
 sig[i,1] <- sum(results$padj < 0.05, na.rm = T)
 sig[i,2] <- sum(results$padj < 0.05 & results$log2FoldChange > 0, na.rm = T)
 sig[i,3] <- sum(results$padj < 0.05 & results$log2FoldChange < 0, na.rm = T)

}

# import number of DEIs from original DEA in up- and down-responders

up <- 6843
down <- 2020

# calculate mean and standard deviation for 10 iterations

mu <- mean(sig[1:10,1])
sdev <- sd(sig[1:10,1])

# calculate standard error

sterr <- sqrt((sdev^2)/10)

# calculate t-statistic

tstat <- (mu - 2016) / sterr

# one-side t-test p-value

1 - pt(tstat,9)

# calculate 97.5% quantile

quantile <- qt(0.975,9)

# calculate upper and lower confidence intervals

CI_lower <- mu - sterr*quantile
CI_upper <- mu + sterr*quantile

# plot results

p_inhouse <- ggplot(data = sig, aes(x=rownames(sig),y=Significant)) +
        geom_col(fill = "#190F9733", color = "#190F97") +
        ggtitle("In-house dataset") +
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
        geom_hline(yintercept = CI_lower, color = "red", size = 0.5, line) +
        geom_hline(yintercept = CI_upper, color = "red", size = 0.5, line) +
        ggpubr::rremove("grid")

p_inhouse





   