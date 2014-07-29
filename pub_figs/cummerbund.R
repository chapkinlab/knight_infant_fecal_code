library(cummeRbund)

datadir <- "/home/jason/repos/sequencing-pipeline/analysis/baby-all-13/cuffdiff"
cuff <- readCufflinks(datadir)

#v<-csVolcano(genes(cuff),"premature","mature", alpha=0.39, showSignificant=TRUE)

s<-csScatter(genes(cuff), "premature","mature", smooth=TRUE) + 
  labs(x="Preterm", y="Term", title=NULL) + 
  theme(text = element_text(size=20, family="sans"))
ggsave("scatter.pdf", plot=s)
#system("convert -density 100 scatter.pdf scatter.png")

#sm <- csScatterMatrix(genes(cuff), replicates=TRUE)
  #labs(x="Preterm", y="Term", title=NULL) + 
  #theme(text = element_text(size=20, family="serif"))
#ggsave("scatter_matrix.png", plot=sm)
#system("convert -density 100 scatter_matrix.pdf scatter_matrix.png")


#myDistHeat<-csDistHeat(genes(cuff),replicates=TRUE)
#print(myDistHeat)

#csd <- csDensity(genes(cuff), replicates=TRUE)
#csb <- csBoxplot(genes(cuff), replicates=TRUE)
