import pylab as p
import pandas as pa
import numpy as np
from scipy import stats as st

p.close('all')

cts = pa.read_excel('11genes.xls', 'ct')
#cts[cts>37] = np.nan
norm = cts.sub(cts.quantile(0.5, axis=1), axis=0)

fid = pa.HDFStore('../rawdata/baby-all-13.h5')
counts = fid['counts']
del counts['type']
logcounts = np.log2(counts)
logcountsnormfac = np.log2(counts.quantile(0.95, axis=1))
logcountsnorm = logcounts.sub(logcountsnormfac, axis=0)

subcounts = logcountsnorm[cts.columns[:-1]].iloc[1:7]
subcounts.values[np.isneginf(subcounts)] = np.nan
subcts = cts.iloc[:,:-1]
subnorm = -norm.iloc[:,:-1]

x,y = (subcounts, subnorm)
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

c1,c2 = 'm','c'
m1,m2 = 'o','D'
#c1,c2 = 'k','k'
#m1,m2 = 'o','D'
#
for i,name in enumerate(pair[0].columns):
    x0,x1 = pair[0].iloc[:3,i], pair[0].iloc[3:,i]
    y0,y1 = pair[1].iloc[:3,i], pair[1].iloc[3:,i]
    #x = pa.Series(x).fillna(-4)
    pt = ax1.plot(getx(i), x0, m1, markersize=ms, color=c1, alpha=0.8)
    ft = ax1.plot(getx(i), x1, m2, markersize=ms, color=c2, alpha=0.8)
    #ft = ax1.plot(getx(i), x1, m2, markersize=ms, color=c2, alpha=0.8, mfc='None')

    ax2.plot(getx(i), y0, m1, markersize=ms, color=c1, alpha=0.8)
    ax2.plot(getx(i), y1, m2, markersize=ms, color=c2, alpha=0.8)
    #ax2.plot(getx(i), y1, m2, markersize=ms, color=c2, alpha=0.8, mfc='None')

    pcrfold = x0.mean() - x1.mean()
    seqfold = y0.mean() - y1.mean()
    print("{}: pcrfold: {:.3f}, seqfold: {:.3f}, cuffdiff: {:.3f}".format(
        name, pcrfold, seqfold, 0.0))
    pcrfoldlist.append(pcrfold)
    seqfoldlist.append(seqfold)

p.legend((pt[0], ft[0]), ('Preterm', 'Term'), loc='best', numpoints=1)
ax1.set_ylabel('RNA-Seq log2(normalized expression)')
ax2.set_ylabel('qPCR log2(Normalized expression)')

ax1.set_xticks(np.arange(pair[0].shape[1]), minor=False)
ax1.set_xticklabels(pair[0].columns, minor=False, rotation=90)
ax2.set_xticks(np.arange(pair[0].shape[1]), minor=False)
ax2.set_xticklabels(()*len(pair[0].columns))
ax1.set_xlim(-0.6, pair[0].shape[1]+0.1)
ax2.set_xlim(-0.6, pair[0].shape[1]+0.1)
ax1.grid(True)
ax2.grid(True)
#ax1.xticks(rotation=90)

#print("averaged: pearson: {:.3f}, spearman: {:.3f}".format(np.mean(pearslist), np.mean(spearlist)))

p.show()
