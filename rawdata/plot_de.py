import os, sys
import pandas as pa
import numpy as np
import matplotlib as mpl
font = {'size' : 8}
mpl.rc('font', **font)
import pylab as p

p.close('all')

if len(sys.argv) != 2:
    print("ERROR: Need h5 file as input argument")
    sys.exit(-1)

fname = sys.argv[1]
fid = pa.HDFStore(fname)

data = fid['counts']
#data = fid['fpkm']
#
del data['type']
datat = data.T[['192_2','196_2','204_2','BF10c','BF1cD','FF12c']]
datat.columns = ['Preterm 1', 'Preterm 2','Preterm 3',
        'Full term 1', 'Full term 2', 'Full term 3']

#thresh = 10
#num = len(datat.columns)
#mat = np.zeros((num,num))
#logicals = datat>thresh
#counts = logicals.sum(axis=0)

#datat = datat[datat.sum(axis=1)>0]
#corrs = datat.corr(method='pearson')
##corrs = datat.corr(method='spearman')
#heatmap(corrs, datat.columns)#, vmin=0)
#p.show()

cdat = pa.read_csv("../cuffdiff/gene_exp.diff", sep='\t', index_col=0)
#qsig = cdat[cdat['significant']=='yes'].sort('p_value')
psig = cdat[cdat['p_value']<0.05].sort('p_value')
#psig_orig = cdat[cdat['p_value']<0.05].sort('p_value')
#craw = pa.read_csv("genes.read_group_tracking", sep='\t', index_col=0)
#psig = psig.sort('log2(fold_change)')

hk_genes = set(map(str.strip, open('HK_only_genes.txt').readlines()))
hk_int_genes = datat.index.intersection(hk_genes)
datat_sub = datat.loc[hk_int_genes]
##Try mean normalizing on the hk genes
datat = datat/datat_sub.mean(axis=0)

# Try median normalizing on the hk genes
# too many zeros to do this naively so we only consider non-zero
#new_hk = datat_sub[(datat_sub>0).prod(axis=1)==1]
#datat = datat/new_hk.median(axis=0)

def distinct(ds1, ds2):
    # Make sure the intervals don't overlap
    h1 = ds1.max(axis=1)
    l1 = ds1.min(axis=1)

    h2 = ds2.max(axis=1)
    l2 = ds2.min(axis=1)
    return ((l1 > h2) + (l2 > h1))

# Now only plot those genes that don't 'overlap'
psig = psig[distinct(datat.ix[psig.index,:3], datat.ix[psig.index,3:])]
N = psig.shape[0]

syms = list('o>')
colors = list('rb')

f = p.figure(figsize=(11,4))
#for i,geneid in enumerate(psig['gene']):
    #for j, samps in enumerate(([0,1,2],[3,4,5])):
        #p.semilogy(i*np.ones(3), datat.ix[geneid,samps]+1e-2, syms[j], color=colors[j], markersize=5, alpha=0.7)

craw = pa.read_csv("../cuffdiff/genes.read_group_tracking", sep='\t', index_col=0)
for i,gene in enumerate(psig['gene']):
    gs = craw.loc[gene,:].groupby('condition')
    for j,lab in enumerate(['mature','premature']):
        g = gs.get_group(lab)
        p.semilogy(i*np.ones(g.shape[0]), g['FPKM']+1e-2, syms[j], color=colors[j], markersize=5, alpha=0.7)

ax = p.gca()
ax.set_xticks(np.arange(N), minor=False)
ax.set_xticklabels(psig.gene_id[:N], minor=False)
p.xticks(rotation=90)
p.grid(True)

#p.ylabel('Normalized counts')
#figname = 'deplot_norm_count_mean.pdf'

p.ylabel('FPKM')
figname = 'deplot_cd_fpkm.pdf'

p.legend(('Preterm','Full term')) 
p.title(figname)
p.savefig(figname, bbox_inches='tight')
