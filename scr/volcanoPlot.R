library(tidyverse)
library(ggrepel)
##### Volcano plot
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <resultsFile> <geneListFile> <experimentLabel> <outputPng> ', call.=FALSE)
}

myd<-read.table(args[1], sep = "\t", header = T)
genes<-read.table(args[2], sep = "\t", header = F)
colnames(genes)<- c("ensemblId", "SYMBOL")
tot<- merge(genes, myd, by = "ensemblId", all =T)
mycol <- c("#171717","#DDDDDD")
tot %>% filter(FDR < 0.5 & PValue < 0.05) %>% mutate(state = ifelse(logFC > -1 & logFC < 1 , "unchanged", "deregulated" )) %>% mutate(trueS = ifelse(is.na(SYMBOL), paste(ensemblId), paste(SYMBOL))) %>% mutate(labels = ifelse(state == "deregulated", paste(trueS), "")) %>% ggplot(aes(logFC, -log10(PValue), color = state)) + geom_point() + theme_bw() + geom_vline(xintercept=c(-1, 1), col="grey") + geom_hline(yintercept=-log10(0.05), col="grey") + scale_color_manual(values = mycol)  + geom_text_repel(aes(label = labels), max.overlaps = 25) + ggtitle(paste(args[3], "FDR <0.5, PValue <0.05, logFC< -1 & > 1", sep = "/"))
ggsave(args[4], heigh = 9, width = 14)