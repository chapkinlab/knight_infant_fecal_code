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

pair = (subcounts, subnorm)
flatpair = map(lambda x: np.array(x).flatten(), pair)

colors = 'rgbmcy'
for i,name in enumerate(pair[0].index):
    #if name != '204_2':
        #continue
    x, y = pair[0].iloc[i], pair[1].iloc[i]
    x = pa.Series(x).fillna(-4)
    p.plot(x, y, '.', markersize=15, label=name, color=colors[i])
    slope, intercept, rval, pval, stderr = st.linregress(x,y)
    p.plot([-4,4], [-4*slope+intercept, slope*4+intercept], '--', color=colors[i])
    print("{}: pearson: {:.3f}, spearman: {:.3f}, slope: {}, intercept: {}".format(
        name, st.pearsonr(x,y)[0], st.spearmanr(x,y)[0],
        slope,intercept))

c = 7
p.plot([-c,c],[-c,c], 'k--')
p.axis([-c,c,-c,c])
x,y = flatpair
print("overall: pearson: {:.3f}, spearman: {:.3f}".format(st.pearsonr(x,y)[0], st.spearmanr(x,y)[0]))

p.xlabel('RNA-Seq log2(count)')
p.ylabel('qPCR log2(expression fold change over median)')
p.legend(loc='best')
p.grid(True)
p.show()


#192_2: pearson: 0.702, spearman: 0.753
#196_2: pearson: 0.108, spearman: 0.546
#204_2: pearson: 0.278, spearman: 0.736
#BF10c: pearson: 0.939, spearman: 0.542
#BF1cD: pearson: 0.670, spearman: 0.309
#FF12c: pearson: 0.619, spearman: 0.627
#overall: pearson: 0.348, spearman: 0.592


