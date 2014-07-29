import pylab as p
import pandas as pa
import numpy as np
from scipy import stats as st

p.close('all')

cts = pa.read_excel('11genes.xls', 'ct')
#cts[cts>37] = np.nan
exp = 2**(40-cts)

#norm = cts.sub(cts['18S'], axis=0)
norm = cts.sub(cts.quantile(0.75, axis=1), axis=0)

fid = pa.HDFStore('../rawdata/baby-all-13.h5')
counts = fid['counts']
#counts = counts.div(counts.median(axis=1))

subcounts = counts[cts.columns[:-1]].iloc[1:7]
subcts = cts.iloc[:,:-1]
subexp = exp.iloc[:,:-1]
subnorm = 2**(-norm.iloc[:,:-1])

pair = (subcounts, subexp)
#pair = (subcounts, subnorm)
flatpair = map(lambda x: np.array(x).flatten(), pair)

#shuffled_index = np.random.permutation(subexp.shape[0])
#np.random.shuffle(flatpair[1])

pearslist, spearlist = [], []

fig = p.figure()
c = 8
colors = 'rgbmcy'
for i,name in enumerate(pair[0].index):
    fig.add_subplot(2,3,i+1)
    #if name != '192_2':
        #continue
    x, y = pair[0].iloc[i], pair[1].iloc[i]
    p.loglog(x+0.1, y+0.1, '.', markersize=15, label=name, color=colors[i])
    #p.plot(x, y, '.', markersize=15, label=name, color=colors[i])
    #slope, intercept, rval, pval, stderr = st.linregress(x,y)
    slope, intercept, rval, pval, stderr = st.linregress(np.log(x+0.1),np.log(0.1+y))
    #print(np.log(x+0.1), np.log(y+0.1))
    xplot = np.logspace(0,4)
    p.loglog(xplot, np.exp(intercept)*xplot**slope, '--', color=colors[i])
    #p.loglog([10,10**3], [10*slope+intercept, slope*10**3+intercept], '--', color=colors[i])
    #p.plot([10,10**3], [10*slope+intercept, slope*10**3+intercept], '--', color=colors[i])
    pears,spear = st.pearsonr(x,y)[0], st.spearmanr(x,y)[0]
    print("{}: pearson: {:.3f}, spearman: {:.3f}, slope: {:.3f}, intercept: {:.3f}, r^2: {:.3f}, pval: {:.3f}".format(
        name, pears, spear, 
        slope,intercept, rval**2, pval))
    pearslist.append(pears)
    spearlist.append(spear)
    p.legend(loc='upper left', numpoints=1)
    p.grid(True)
    #p.plot([-c,c],[-c,c], 'k--')
    #p.axis([-4,5,-5,8])

x,y = flatpair
print("overall: pearson: {:.3f}, spearman: {:.3f}".format(st.pearsonr(x,y)[0], st.spearmanr(x,y)[0]))
print("averaged: pearson: {:.3f}, spearman: {:.3f}".format(np.mean(pearslist), np.mean(spearlist)))

fig.text(0.4, 0.04, 'RNA-Seq counts', fontsize=14)
fig.text(0.06, 0.8, 'qPCR un-normalized expression 2^(40-ct)', rotation='vertical', fontsize=14)
p.show()

#192_2: pearson: 0.702, spearman: 0.753
#196_2: pearson: 0.108, spearman: 0.546
#204_2: pearson: 0.278, spearman: 0.736
#BF10c: pearson: 0.939, spearman: 0.542
#BF1cD: pearson: 0.670, spearman: 0.309
#FF12c: pearson: 0.619, spearman: 0.627
#overall: pearson: 0.348, spearman: 0.592


