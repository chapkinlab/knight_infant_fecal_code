library(reshape2)
library(ggplot2)
library(ShortRead)
library(scales)

# Read in ERCC counts for each sample
erccdat <- read.csv('../rawdata/baby-all-13-ercc.csv', header=T, row.names=1)

# Cleanup
erccdat[is.na(erccdat)] <- 0
colnames(erccdat) <- sapply(colnames(erccdat), function(x) gsub(".","-",x,fixed=T))
types <- erccdat$type
levels(types)[levels(types)=="premature"] <- "Preterm"
levels(types)[levels(types)=="fullterm"] <- "Term"
erccdat$type <- NULL
erccdat <- as.data.frame(t(erccdat))
erccdat["192_1"] <- NULL
erccdat$s_8_sequence <- NULL
erccdat$s_7_sequence <- NULL
colnames(erccdat) <- c("Preterm 1","Preterm 2","Preterm 3", "Term 1","Term 2","Term 3")#,"Full Pooled")

ercc.cors <- cor(erccdat, method='pearson')
ercc.scors <- cor(erccdat, method='spearman')

# Read in ERCC Lengths
erccdat$lengths <- width(readFasta("/mnt/datab/refs/igenomes/ercc/ERCC92.fa"))

# Read in ERCC concentrations
dat <- read.csv("/mnt/datab/refs/igenomes/ercc/cms_095046.txt", sep="\t", row.names=2)[3]
names(dat) <- "conc"
erccdat <- merge(erccdat, dat, by=0)
rownames(erccdat) <- erccdat$Row.names
erccdat$Row.names <- NULL

# Melt
vtall <- melt(erccdat, id.vars=c("lengths","conc"), value.name="reads", variable.name="sample")

# Add back in types
vtall$type <- types[match(vtall$sample, colnames(erccdat))]

# Normalize for length
vtall$readsnorm <- vtall$reads/vtall$lengths

vtall$readsnorm[vtall$readsnorm==0] <- 0.1
vtall$reads[vtall$reads==0] <- 0.1

# Calculate molecules
vtall$mols <- vtall$conc/100*5*6.022e23*1e-18

fig <- ggplot(vtall, aes(x=mols, y=reads, color=sample)) + 
  #facet_wrap( ~ sample ) +
  stat_smooth(method='lm', se=F) +
  geom_point(aes(x=mols, y=reads, shape=type, size=lengths, color=sample), alpha=0.8) +
  scale_x_log10(breaks=10^(0:9), labels=trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks=10^(-2:8), labels=trans_format("log10", math_format(10^.x))) +
  xlab("Number of molecules") + 
  ylab("Observed reads") 
  #theme_bw()
#scale_x_log10() + scale_y_log10()
print(fig)
ggsave("ERCC_mols.pdf", plot=fig)

