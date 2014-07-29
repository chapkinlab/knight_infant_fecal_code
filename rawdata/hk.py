import pandas as pa
import pylab as p

fid = pa.HDFStore('baby-all-13.h5')
dat = fid['counts']

hk_genes = set(map(str.strip, open('HK_only_genes.txt').readlines()))

hk_int_genes = hk_genes.intersection(dat.columns)

dat_sub = dat.loc[:,hk_int_genes]
del dat['type']

# Try normalizing on the hk genes
dat = (dat.T/dat_sub.mean(axis=1)).T

# Nahh.. doesn't seem to help...
# woah, these pearson correlations are interesting....
#corrs = dat.T.corr(method='spearman')
corrs = dat.T.corr(method='pearson')
p.pcolor(corrs, cmap=p.cm.RdYlBu, vmin=-1.0, vmax=1.0)
p.colorbar()
p.show()
