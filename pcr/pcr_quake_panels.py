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
subnorm.index = ['Preterm 1','Preterm 2','Preterm 3','Term 1','Term 2','Term 3']

pair = (subcounts, subnorm)
flatpair = map(lambda x: np.array(x).flatten(), pair)

pearslist, spearlist, slopelist, rlist = [], [], [], []

fig = p.figure()
c = 8
colors = 'rgbmcy'
for i,name in enumerate(pair[1].index):
    ax = fig.add_subplot(2,3,i+1)
    #if name != '204_2':
        #continue
    x, y = pair[0].iloc[i], pair[1].iloc[i]
    x = pa.Series(x).fillna(-4)
    ax.plot(x, y, '.', markersize=15, label=name, color=colors[i])
    slope, intercept, rval, pval, stderr = st.linregress(x,y)
    p.plot([-5,8], [-5*slope+intercept, slope*8+intercept], '--', color=colors[i])
    pears,spear = st.pearsonr(x,y)[0], st.spearmanr(x,y)[0]
    print("{}: pearson: {:.3f}, spearman: {:.3f}, slope: {:.3f}, intercept: {:.3f}, r^2: {:.3f}, pval: {:.3f}".format(
        name, pears, spear, 
        slope,intercept, rval**2, pval))
    ax.text(0.1,0.7, "pearson: {:.3f}\nspearman: {:.3f}\nslope: {:.3f}\nintercept: {:.3f}\nr^2: {:.3f}\npval: {:.3f}".format(
        pears, spear, 
        slope,intercept, rval**2, pval), fontsize=12,
        #horizontalalignment='center',
         verticalalignment='center',
         transform=ax.transAxes)
    pearslist.append(pears)
    spearlist.append(spear)
    slopelist.append(slope)
    rlist.append(rval**2)
    p.legend(loc='upper left', numpoints=1, fontsize=12)
    #p.grid(True)
    p.plot([-c,c],[-c,c], 'k--')
    p.axis([-4,5,-5,8])

x,y = flatpair
x = pa.Series(x).fillna(-4)
print("overall: pearson: {:.3f}, spearman: {:.3f}".format(st.pearsonr(x,y)[0], st.spearmanr(x,y)[0]))
print("averaged: pearson: {:.3f}, spearman: {:.3f}, slope: {:.3f}, r^2: {:.3f}".format(*map(np.mean, [pearslist, spearlist, slopelist,rlist])))


fig.text(0.4, 0.04, 'RNA-Seq log2(count)', fontsize=14)
fig.text(0.06, 0.8, 'qPCR log2(expression fold change over median)', rotation='vertical', fontsize=14)
p.show()


#192_2: pearson: 0.702, spearman: 0.753
#196_2: pearson: 0.108, spearman: 0.546
#204_2: pearson: 0.278, spearman: 0.736
#BF10c: pearson: 0.939, spearman: 0.542
#BF1cD: pearson: 0.670, spearman: 0.309
#FF12c: pearson: 0.619, spearman: 0.627
#overall: pearson: 0.348, spearman: 0.592


