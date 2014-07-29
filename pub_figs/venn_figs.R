library(reshape2)
library(ggplot2)
library(VennDiagram)

# Read in mapped counts for each sample
rawcts <- read.csv('../rawdata/baby-all-13.fpkm.T.csv', header=T, row.names=1, stringsAsFactors=F)
key <- read.csv('../rawdata/baby-key.csv', header=T)

x <- lapply(as.data.frame(rawcts[c("BF10c","BF1cD","FF12c","s_7_sequence")]>10), which)
names(x) <- c("Full term 1", "Full term 2", "Full term 3", "Pooled full term")
fig <- venn.diagram(x, fill=c("orange","green","blue","red"), 
                    alpha=c(0.3,0.3,0.3,0.3), filename = "fullterm_rpkm_10_venn.pdf", 
                    fontfamily=rep("sans",15), cat.fontfamily=rep("sans",4))

system("convert -density 90 fullterm_rpkm_10_venn.pdf fullterm_rpkm_10_venn.pdf.png")
