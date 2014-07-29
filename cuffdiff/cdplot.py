import os, sys
import pandas as pa
import pylab as p
import numpy as np
from scipy import stats
#from statsmodels.graphics.boxplots import violinplot
#import rpy2.robjects as ro
#r = ro.r

p.close('all')

if len(sys.argv) != 2:
    print("ERROR: Need h5 file as input argument")
    sys.exit(-1)

fname = sys.argv[1]

fid = pa.HDFStore(fname)

data = fid['fpkm']
datat = data.T

cdat = pa.read_csv("gene_exp.diff", sep='\t', index_col=0)

qsig = cdat[cdat['significant']=='yes'].sort('p_value')
psig = cdat[cdat['p_value']<0.05].sort('p_value')
psig_orig = cdat[cdat['p_value']<0.05].sort('p_value')

#sys.exit()


craw = pa.read_csv("genes.read_group_tracking", sep='\t', index_col=0)

syms = list('>o')
colors = list('br')

N = psig.shape[0]
#psig = psig.iloc[:N,:]
psig = psig.sort('log2(fold_change)')

#inds = np.hstack(( np.arange(0,N/2), np.arange(psig.shape[0]-N/2, psig.shape[0]) ))
#psig = psig.iloc[inds,:]
f = p.figure(figsize=(70,8))

for i,gene in enumerate(psig.gene_id[:N]):
    gs = craw.loc[gene,:].groupby('condition')
    for j,lab in enumerate(['premature','mature']):
        g = gs.get_group(lab)
        p.semilogy(i*np.ones(g.shape[0]), g['FPKM']+1e-2, syms[j], color=colors[j], markersize=10, alpha=0.7)

ax = p.gca()
ax.set_xticks(np.arange(N), minor=False)
ax.set_xticklabels(psig.gene_id[:N], minor=False)
p.xticks(rotation=75)
p.ylabel('FPKM')
p.grid(True)
p.legend(('Full term','Preterm')) # NB: This inversion is intentional
#p.savefig('deplot_cd_fpkm.pdf', bbox_inches='tight')
p.show()
    
