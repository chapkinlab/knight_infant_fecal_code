import os, sys
import subprocess as sb
import pandas as pa
import numpy as np
from scipy import stats
import statsmodels.graphics.boxplots
#import rpy2.robjects as ro
#r = ro.r
import matplotlib as mpl
font = {'size' : 18}
mpl.rc('font', **font)
import pylab as p

p.close('all')

if len(sys.argv) != 2:
    print("ERROR: Need h5 file as input argument")
    sys.exit(-1)

fname = sys.argv[1]

def heatmap(corrs, labels, vmin=-1, vmax=1):
    fig, ax = p.subplots(figsize=(10,8))
    ax.set_frame_on(False)
    #p.pcolor(corrs, cmap=p.cm.OrRd, vmin=vmin, vmax=vmax)
    #p.pcolor(corrs, cmap=p.cm.RdYlBu, vmin=vmin, vmax=vmax)
    #p.pcolor(corrs, cmap=p.cm.Greys, vmin=vmin, vmax=vmax)
    p.imshow(corrs, cmap=p.cm.RdYlBu, vmin=vmin, vmax=vmax, interpolation='none')
    #p.imshow(corrs, cmap=p.cm.Greys, vmin=vmin, vmax=vmax, interpolation='none')
    #p.imshow(corrs, cmap=p.cm.Greys, vmin=vmin, vmax=vmax, interpolation='none')

    ## for pcolor:
    ## put the major ticks at the middle of each cell
    #ax.set_xticks(np.arange(corrs.shape[0])+0.5, minor=False)
    #ax.set_yticks(np.arange(corrs.shape[1])+0.5, minor=False)

    ## want a more natural, table-like display
    #ax.invert_yaxis()
    #ax.xaxis.tick_top()

    # for imshow:
    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(corrs.shape[0]), minor=False)
    ax.set_yticks(np.arange(corrs.shape[1]), minor=False)

    # want a more natural, table-like display
    #ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(datat.columns, minor=False)
    ax.set_yticklabels(datat.columns, minor=False)

    p.xticks(rotation=90)

    ax.grid(False)
    # Turn off all the ticks
    ax = p.gca()

    for t in ax.xaxis.get_major_ticks(): 
        t.tick1On = False 
        t.tick2On = False 
    for t in ax.yaxis.get_major_ticks(): 
        t.tick1On = False 
        t.tick2On = False  

    p.colorbar()

fid = pa.HDFStore(fname)
colors = list('bggbmyckr')

data = fid['counts']
#data = fid['fpkm']
del data['type']
datat = data.T[['192_2','196_2','204_2','BF10c','BF1cD','FF12c', 's_7_sequence']]
datat.columns = ['Preterm 1', 'Preterm 2','Preterm 3',
        'Term 1', 'Term 2', 'Term 3', 'Pooled term']

hk_genes = set(map(str.strip, open('HK_only_genes.txt').readlines()))
hk_int_genes = datat.index.intersection(hk_genes)
datat_sub = datat.loc[hk_int_genes]
##Try mean normalizing on the hk genes
datat = datat/datat_sub.mean(axis=0)

## Try median normalizing on the hk genes
## too many zeros to do this naively so we only consider non-zero
##new_hk = datat_sub[(datat_sub>0).prod(axis=1)==1]
##datat = datat/new_hk.median(axis=0)

##thresh = 10
##num = len(datat.columns)
##mat = np.zeros((num,num))
##logicals = datat>thresh
##counts = logicals.sum(axis=0)

##for i in range(num):
    ##for j in range(num):
        ##if i == j:
            ##temp = counts[i]*1.0 / counts[i]
            ###temp = counts[i]
        ##else:
            ##temp = (logicals.iloc[:,i] * logicals.iloc[:,j]).sum()*1.0 / counts[j]
            ###temp = (logicals.iloc[:,i] * logicals.iloc[:,j]).sum()
        ##mat[i,j] = temp

##heatmap(mat, datat.columns, vmin=0, vmax=mat.max())
##for i in range(num):
    ##for j in range(num):
        ##if i == j:
            ###p.text(i+0.1,j+0.6, "%d"%counts[i])
            ##p.text(i+0.3,j+0.55, "%d"%counts[i], color='w')
        ##else:
            ###p.text(i+0.1,j+0.6, "%d"% mat[j,i])
            ##p.text(i+0.4,j+0.55, "%d"% int(100*mat[j,i]))
            ###p.text(i+0.2,j+0.7, "%.2f"%mat[i,j])

##y = 9.8
##y = 6.5
##p.text(0,y, """Number at i,j: 'i has __ percentage of j's expressed genes'
##when considering a threshold >{} FPKM""".format(thresh))
##p.show()
##sys.exit()

#datat = datat[datat.sum(axis=1)>0]
#corrs = datat.corr(method='pearson')
#corrs = datat.corr(method='spearman')
#heatmap(corrs, datat.columns)#, vmin=0)
#p.tight_layout()
#figname = 'spearman.pdf'
#p.savefig(figname, bbox_inches='tight')
#sb.check_call('convert -density 100 {} {}.png'.format(figname,figname), shell=True)
#p.show()
#sys.exit()

#data = fid['counts']
#datat = data.T
######################################################
######################################################
datat[np.log10(datat).min(axis=1)<-30] = -30

statsmodels.graphics.boxplots.violinplot([np.log10(datat[c][datat[c]>0]) for 
    c in datat.columns], labels=datat.columns, plot_opts={'label_rotation':30})
p.grid(True)
p.ylabel('log10(counts)')
p.title('Expression Histograms - Counts')
p.show()

######################################################
######################################################
p.figure()
for i,c in enumerate(datat.columns):
    #kde = stats.gaussian_kde(np.log(datat[c][datat[c]>0]))
    #kde = stats.gaussian_kde(datat[c][datat[c]>0])
    #kde = stats.gaussian_kde(datat[c])
    #x = np.linspace(0,3,100)
    #x = np.logspace(0,3,100)
    #p.plot(x,kde(x), color=colors[i], label=c)
    #p.semilogx(x,kde(x), color=colors[i], label=c)
    bins = np.logspace(0,5,20)
    cs,_ = np.histogram(datat[c],bins=bins)
    p.loglog(bins[:-1],cs, color=colors[i], label=c)
    
p.legend(fontsize=8, loc='best')
p.ylabel('# of genes')
p.xlabel('counts')
p.grid(True)
p.title('Smoothed Expression Histograms - Counts')
######################################################
######################################################
#p.figure()
#[datat[c].hist(log=True, bins=np.logspace(0,5,20), histtype='step', label=c, 
    #color=colors[i]) for (i,c) in enumerate(datat.columns)]
#p.xscale('log')
#p.legend(fontsize=8)
#p.ylabel('# of genes')
#p.xlabel('counts')

######################################################
######################################################
data = fid['fpkm']
datat = data.T
statsmodels.graphics.boxplots.violinplot([np.log10(datat[c][datat[c]>1e-3]) for 
    c in datat.columns], labels=datat.columns, plot_opts={'label_rotation':30})
p.grid(True)
p.ylabel('log10(fpkm)')
p.title('Expression Histograms - FPKM')

######################################################
######################################################
p.figure()
for i,c in enumerate(datat.columns):
    filt = np.log(datat[c][datat[c]>0])
    ##filt = np.log(datat[c].fillna(1e-30))
    #kde = stats.gaussian_kde(filt)
    #x = np.linspace(-8,8,100)
    #p.plot(x,kde(x), color=colors[i], label=c)
    bins = np.linspace(-3,5,30)
    cs,_ = np.histogram(filt,bins=bins)
    p.semilogy(bins[:-1],cs, color=colors[i], label=c)
p.legend(fontsize=8, loc='best')
p.xlabel('log(fpkm)')
p.ylabel('# of genes')
p.grid(True)
p.title('Smoothed Expression Histograms - FPKM')

######################################################
######################################################
#p.figure()
#[datat[c].hist(log=True, bins=np.logspace(-3,4,20), histtype='step', label=c,
#color=colors[i]) for (i,c) in enumerate(datat.columns)]
#p.xscale('log')
#p.legend(fontsize=8)
#p.xlabel('fpkm')
#p.ylabel('# of genes')

p.show()
