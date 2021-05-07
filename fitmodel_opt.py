# This code fits the gmchiicorr model to a measured small-scale 2-point cross correlation function

import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
import time
from scipy.optimize import curve_fit
#from scipy.optimize import curve_fit
#import emcee
#import pickle
#from multiprocessing import Pool
#from multiprocessing import cpu_count
#import os
#os.environ["OMP_NUM_THREADS"] = "1"
#import corner

from gmchiicorr import eval_w

#np.random.seed(666)

# ## Name of galaxy

galaxy=sys.argv[1]
print(galaxy)


# ## Read GMC and Random catalogs coordinates

xygmc=ascii.read('./output/'+galaxy+'_xy_gmc.txt')
x1=xygmc['col1'].data
y1=xygmc['col2'].data

xyrand=ascii.read('./output/'+galaxy+'_xy_rand.txt')
xr=xyrand['col1'].data
yr=xyrand['col2'].data

# ## Read Observed Correlation Function

obscorr=ascii.read('./output/'+galaxy+'_corr_small.txt')
r0obs=obscorr['col1'].data
w0obs=obscorr['col2'].data
ew0obs=obscorr['col3'].data

## Read fraction of HII regions to GMCs

fhgtab=ascii.read('./output/'+galaxy+'_fhg.txt')
fhgobs=fhgtab['col2'].data[0]   
efhgobs=fhgtab['col2'].data[1]

# # Define Priors and Likelihood Functions

# Select which range in r to fit
selr=(r0obs<=500)

# First fit correlation function and then fhg
#print("Fitting Correlation Function")
#def func1(r, p0, p1, p2, p3, p4, p5, p6):
#    bins=r[0:-1]
#    r0, w0, ew0, fhg0 = eval_w(l0=p0, rc0=p1, tc0=p2, ts0=p3, tfb0=p4, Ng0=p5, voff0=p6, bins=bins, Nsamples=150)  #Nsamples=150 yields rms smaller than measurement errors
#    return np.concatenate([w0,np.array([fhg0])])    
#p0=np.array([200, 50, 20, 10, 2, 5, 5])
#auxr=np.concatenate([r0obs,np.array([-1])])
#auxw=np.concatenate([w0obs,np.array([fhgobs])])
#auxew=np.concatenate([ew0obs,np.array([efhgobs])])
#pstart, pcov = curve_fit(func1, auxr, auxw, p0=p0, sigma=auxew, method='lm', epsfcn=0.01)
#pbest=pstart

# Fixing ts=5
def func1(r, p0, p1, p2, p4, p5, p6):
    bins=r[0:-1]
    r0, w0, ew0, fhg0 = eval_w(l0=p0, rc0=p1, tc0=p2, ts0=5, tfb0=p4, Ng0=p5, voff0=p6, bins=bins, Nsamples=150)  #Nsamples=150 yields rms smaller than measurement errors
    return np.concatenate([w0,np.array([fhg0])])    
p0=np.array([200, 50, 10, 2, 5, 10])
auxr=np.concatenate([r0obs,np.array([-1])])
auxw=np.concatenate([w0obs,np.array([fhgobs])])
auxew=np.concatenate([ew0obs,np.array([efhgobs])])
pstart, pcov = curve_fit(func1, auxr, auxw, p0=p0, sigma=auxew, method='lm', epsfcn=0.01)
pstart=np.array([pstart[0], pstart[1], pstart[2], 5.0, pstart[3], pstart[4], pstart[5]])
pbest=pstart


print("Best-fit Parameters:", pbest)
r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], bins=r0obs)

## Make best fit plot

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(r0, w0, '-o', color='green', alpha=1.0, label="fhg="+"{:.2f}".format(fhg0))    
ax.errorbar(r0obs, w0obs, ew0obs, fmt="o", color='black', capsize=5, alpha=0.5)
ax.plot(r0obs, w0obs, 'o', color='black', alpha=0.5)
#ax.set_xlim(0, 00)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.set_title(galaxy+" ; fhg="+"{:.2f}".format(fhgobs)+" ("+"{:.2f}".format(efhgobs)+")", fontsize=30)
ax.tick_params(labelsize=20)
ax.plot([0], [0], color='white', label="l="+"{:.2f}".format(pbest[0]))
ax.plot([0], [0], color='white', label="rc="+"{:.2f}".format(pbest[1]))
ax.plot([0], [0], color='white', label="tc="+"{:.2f}".format(pbest[2]))
ax.plot([0], [0], color='white', label="ts="+"{:.2f}".format(pbest[3]))
ax.plot([0], [0], color='white', label="tfb="+"{:.2f}".format(pbest[4]))
ax.plot([0], [0], color='white', label="Ng="+"{:.2f}".format(pbest[5]))
ax.plot([0], [0], color='white', label="voff="+"{:.2f}".format(pbest[6]))
ax.legend(fontsize=20)
plt.savefig('./plots/'+galaxy+'_fitopt.png')
