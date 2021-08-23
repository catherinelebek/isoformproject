# This script is to visualise the overlap between BP GO terms enriched among DEIs from the various paired-sample isoform-level DEA analyses

library(VennDiagram) # for plotting venn diagrams
library(ggvenn) # for formatting venn diagrams

up.up <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/up_up.txt")
down.down <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/down_down.txt")
up.up.jarid2.tss <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/up_up_jarid2_tss.txt")
up.up.jarid2.ass <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/up_up_jarid2_ass.txt")
up.up.jarid2.onlyass <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/up_up_jarid2_onlyass.txt")
up.up.jarid2.notass <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/up_up_jarid2_notass.txt")
up.down <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/up_down.txt")
down.down.jarid2.tss <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/down_down_jarid2_tss.txt")
down.down.jarid2.ass <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/down_down_jarid2_ass.txt")
down.down.jarid2.notass <- read.table("~/Documents/Semester3/Project/Results/resultsanalysis/GO_DEA/GO_Terms/down_down_jarid2_notass.txt")

up.up <- c(t(up.up))
down.down <- c(t(down.down))
up.up.jarid2.tss <- c(t(up.up.jarid2.tss))
up.up.jarid2.ass <- c(t(up.up.jarid2.ass))
up.up.jarid2.onlyass <- c(t(up.up.jarid2.onlyass))
up.up.jarid2.notass <- c(t(up.up.jarid2.notass))
up.down <- c(t(up.down))
down.down.jarid2.tss <- c(t(down.down.jarid2.tss))
down.down.jarid2.ass <- c(t(down.down.jarid2.ass))
down.down.jarid2.notass <- c(t(down.down.jarid2.notass))

total <- c(up.up,down.down,up.up.jarid2.tss,up.up.jarid2.ass, up.up.jarid2.onlyass, up.up.jarid2.notass, up.down, down.down.jarid2.tss, down.down.jarid2.ass, down.down.jarid2.notass)
total <- unique(total)

go_res <- data.frame("GO_term" = total)
go_res$up.up <- go_res$GO_term %in% up.up
go_res$down.down <- go_res$GO_term %in% down.down
go_res$up.up.jarid2.tss <- go_res$GO_term %in% up.up.jarid2.tss
go_res$up.up.jarid2.ass <- go_res$GO_term %in% up.up.jarid2.ass
go_res$up.up.jarid2.onlyass <- go_res$GO_term %in% up.up.jarid2.onlyass
go_res$up.up.jarid2.notass <- go_res$GO_term %in% up.up.jarid2.notass
go_res$up.down <- go_res$GO_term %in% up.down
go_res$down.down.jarid2.tss <- go_res$GO_term %in% down.down.jarid2.tss
go_res$down.down.jarid2.ass <- go_res$GO_term %in% down.down.jarid2.ass
go_res$down.down.jarid2.notass <- go_res$GO_term %in% down.down.jarid2.notass

go_res

# plotting overlap in GO terms between up.up and down.down

ggplot(go_res, aes(A = up.up, B = down.down)) +
  geom_venn(set_names = c("Up in URs","Down in DRs")) +
  theme_bw() +
  ggpubr::rremove("grid") +
  ggpubr::rremove("axis")


# plotting overlap in GO terms between up.up and down.down

ggplot(go_res, aes(A = up.up.jarid2.tss, B = up.up.jarid2.notass)) +
  geom_venn(set_names = c("JARID2 TSS","Not JARID2")) +
  theme_bw() +
  ggpubr::rremove("grid") +
  ggpubr::rremove("axis")

ggplot(go_res, aes(A = up.up.jarid2.tss, B = up.up.jarid2.onlyass)) +
  geom_venn(set_names = c("JARID2 TSS","Not JARID2")) +
  theme_bw() +
  ggpubr::rremove("grid") +
  ggpubr::rremove("axis")

ggplot(go_res, aes(A = up.up.jarid2.tss, B = up.up.jarid2.notass)) +
  geom_venn(set_names = c("JARID2 TSS","Not JARID2")) +
  theme_bw() +
  ggpubr::rremove("grid") +
  ggpubr::rremove("axis")


ggplot(go_res, aes(A = down.down.jarid2.tss, B = down.down.jarid2.notass)) +
  geom_venn(set_names = c("JARID2 TSS","Not JARID2")) +
  theme_bw() +
  ggpubr::rremove("grid") +
  ggpubr::rremove("axis")


