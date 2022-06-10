library(tidyverse)
library(ggrepel)
library(gplots)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
 stop('specify five arguments: <resultsFile> <geneListFile> <experimentLabel> <VolcanoPlotPng> <countsMatrixFile> <HeatmapPlot> ', call.=FALSE)
}


myd<-read.table(args[1], sep = "\t", header = T)
genes<-read.table(args[2], sep = "\t", header = F)
colnames(genes)<- c("ensemblId", "SYMBOL")
tot<- merge(genes, myd, by = "ensemblId", all =T)
mycol <- c("#171717","#DDDDDD")
tot %>% filter(FDR < 0.05 & PValue < 0.05) %>% mutate(state = ifelse(logFC > -1.5 & logFC < 1.5 , "unchanged", "deregulated" )) %>% mutate(trueS = ifelse(is.na(SYMBOL), paste(ensemblId), paste(SYMBOL))) %>% mutate(labels = ifelse(state == "deregulated", paste(trueS), "")) %>% ggplot(aes(logFC, -log10(PValue), color = state)) + geom_point() + theme_bw() + geom_vline(xintercept=c(-1.5, 1.5), col="grey") + geom_hline(yintercept=-log10(0), col="grey") + scale_color_manual(values = mycol)  + geom_text_repel(aes(label = labels), max.overlaps = 25) + ggtitle(paste(args[3], "FDR <0.05, PValue <0.05, logFC< -1 & > 1", sep = "/"))
ggsave(args[4], heigh = 9, width = 14)

counts<-read.table(args[5], sep = "\t", header =T)
new<-tot %>% filter(FDR < 0.05 & PValue < 0.05) %>% mutate(state = ifelse(logFC > -1.5 & logFC < 1.5 , "unchanged", "deregulated" )) %>% mutate(trueS = ifelse(is.na(SYMBOL), paste(ensemblId), paste(SYMBOL))) %>% mutate(labels = ifelse(state == "deregulated", paste(trueS), ""))
final<-merge(new, counts, by = "ensemblId")
mydf<-final %>% filter(abs(logFC) > 1.5 | abs(logFC)< -1.5 & FDR<0.05 & PValue < 0.05) %>% select(trueS, 12:17)
mat<-as.matrix(mydf[,2:7])
rownames(mat) <- mydf$trueS
png(args[6], width = 18, height = 30, units = "cm", res = 300)
heatmap.2(mat, dendrogram='none', Rowv=TRUE, Colv=FALSE, reorderfun=function(d,w) reorder(d, w, agglo.FUN=mean), distfun=function(x) as.dist(1-cor(t(x))), hclustfun=function(x) hclust(x, method="complete"),trace='none',col = greenred(100),scale='row',cexRow=0.7, cexCol = 0.8, offsetCol = 0.5, lwid = c(5,15), lhei = c(3,15), margins = c(6, 10))
dev.off()