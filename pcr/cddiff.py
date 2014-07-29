import sys
import pylab as p
import pandas as pa
import numpy as np
from scipy import stats as st

p.close('all')

cts = pa.read_excel('11genes.xls', 'ct')
cts.index = ['Preterm 1','Preterm 2','Preterm 3','Full term 1','Full term 2','Full term 3']
cts = cts.sort_index(axis=0)
#cts[cts>37] = np.nan
norm = cts.sub(cts.quantile(0.5, axis=1), axis=0)


craw = pa.read_csv("../cuffdiff/genes.read_group_tracking", sep='\t', index_col=0)
indname = craw.index.name
craw = craw.reset_index()
cddata = pa.pivot_table(craw, rows=indname, cols=['condition','replicate'], values='FPKM')
cddata = cddata.T.loc[:,norm.columns[:-1]]
cddata = cddata.reset_index()
cddata.loc[cddata['condition']=='mature','condition'] = 'Preterm' ##NB: This is intentional as the results
#from cufflinks are actually flipped
cddata.loc[cddata['condition']=='premature','condition'] = 'Full term'
cddata.condition = cddata.condition + cddata.replicate.map(str)
del cddata['replicate']
cddata = cddata.set_index('condition')
#cddata.index = cddata.index.set_levels([['Full term','Preterm'],[0,1,2]])
#cddata.index = map(lambda x: ' '.join(map(str,x)),cddata.index.to_series())
cddata = np.log2(cddata.sort_index(axis=0))

subnorm = -norm.iloc[:,:-1]

x,y = (cddata, subnorm)
xsort = x.mean(axis=0)
#xsort = x.iloc[:3].mean(axis=0)-x.iloc[3:].mean(axis=0)
xsort.sort(ascending=False)
pair = (x.reindex_axis(xsort.index, axis=1), y.reindex_axis(xsort.index, axis=1))

pcrfoldlist, seqfoldlist = [], []

fig = p.figure()
ax1 = fig.add_subplot(2,1,1)
ax2 = fig.add_subplot(2,1,2)
ms = 10
def getx(i):
    return i+(np.random.rand(3)-0.5)*0.3
    #return i+np.random.rand(3)*0.1

for i,name in enumerate(pair[0].columns):
    x0,x1 = pair[0].iloc[:3,i], pair[0].iloc[3:,i]
    y0,y1 = pair[1].iloc[:3,i], pair[1].iloc[3:,i]
    #x = pa.Series(x).fillna(-4)
    ft = ax1.plot(getx(i), x0, 'o', markersize=ms, color='k', alpha=0.8)
    pt = ax1.plot(getx(i), x1, 'D', markersize=ms, color='k', alpha=0.8, mfc='None')

    ax2.plot(getx(i), y0, 'o', markersize=ms, color='k', alpha=0.8)
    ax2.plot(getx(i), y1, 'D', markersize=ms, color='k', alpha=0.8, mfc='None')

    pcrfold = x0.mean() - x1.mean()
    seqfold = y0.mean() - y1.mean()
    print("{}: pcrfold: {:.3f}, seqfold: {:.3f}".format(
        name, pcrfold, seqfold))
    pcrfoldlist.append(pcrfold)
    seqfoldlist.append(seqfold)

p.legend((pt[0], ft[0]), ('Preterm', 'Full term'), loc='best', numpoints=1)
ax1.grid(True)
ax2.grid(True)
ax1.set_ylabel('RNA-Seq log2(FPKM)')
ax2.set_ylabel('qPCR log2(Normalized expression)')

ax1.set_xticks(np.arange(pair[0].shape[1]), minor=False)
ax1.set_xticklabels(pair[0].columns, minor=False, rotation=90)
ax2.set_xticklabels(())
ax1.set_xlim(-0.5, pair[0].shape[1])
ax2.set_xlim(-0.5, pair[0].shape[1])
#ax1.xticks(rotation=90)

#print("averaged: pearson: {:.3f}, spearman: {:.3f}".format(np.mean(pearslist), np.mean(spearlist)))

p.show()
