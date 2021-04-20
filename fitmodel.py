# This code fits the gmchiicorr model to a measured small-scale 2-point cross correlation function

import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
import time
from scipy.optimize import curve_fit
import emcee
import pickle
from multiprocessing import Pool
from multiprocessing import cpu_count
import os
os.environ["OMP_NUM_THREADS"] = "1"
import corner

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

# Define prior parameter space
lrange=np.array([20,500])
rcrange=np.array([5,150])
tcrange=np.array([1,50])
tsrange=np.array([1,50])
tfbrange=np.array([1,50])
Ngrange=np.array([1,30])
voffrange=np.array([0,50])

# Define uniform prior distribution
def log_prior(p):
    l1, rc1, tc1, ts1, tfb1, Ng1, voff1 = p
    if lrange[0]<l1<lrange[1] and rcrange[0]<rc1<rcrange[1] and tcrange[0]<tc1<tcrange[1] and tsrange[0]<ts1<tsrange[1] and tfbrange[0]<tfb1<tfbrange[1] and Ngrange[0]<Ng1<Ngrange[1] and voffrange[0]<voff1<voffrange[1]:
        return 0.0
    else:
        return -np.inf

# Define likelihood*prior distriibution
def log_prob(p):
    lprior=log_prior(p)
    if np.isfinite(lprior):
        r0, w0, ew0, fhg0 = eval_w(l0=p[0], rc0=p[1], tc0=p[2], ts0=p[3], tfb0=p[4], Ng0=p[5], voff0=p[6], bins=r0obs, Nsamples=150)  #Nsamples=150 yields rms smaller than measurement errors
        res=w0-w0obs[selr]
        sig=ew0obs[selr]
        prob=1/(2*np.pi*sig**2)*np.exp(-0.5*(res/sig)**2)

        resfhg=fhg0-fhgobs
        sigfhg=efhgobs
        probfhg=1/(2*np.pi*sigfhg**2)*np.exp(-0.5*(resfhg/sigfhg)**2)

        logp=lprior+np.sum(np.log(prob))+np.log(probfhg)

        if not np.isfinite(logp):
            return -np.inf
        return logp
    else:
        return -np.inf     


# # Set up MCMC

ndim=7
nwalkers=16

p0 = np.zeros((nwalkers, ndim))
p0[:,0]=np.random.uniform(lrange[0], lrange[1], nwalkers)
p0[:,1]=np.random.uniform(rcrange[0], rcrange[1], nwalkers)
p0[:,2]=np.random.uniform(tcrange[0], tcrange[1], nwalkers)
p0[:,3]=np.random.uniform(tsrange[0], tsrange[1], nwalkers)
p0[:,4]=np.random.uniform(tfbrange[0], tfbrange[1], nwalkers)
p0[:,5]=np.random.uniform(Ngrange[0], Ngrange[1], nwalkers)
p0[:,6]=np.random.uniform(voffrange[0], voffrange[1], nwalkers)

# # Run MCMC Chain

Nmc=2000
#Nmc=2

print("Starting MCMC")
t0full=time.time()

## run MCMC with multiple cpus
with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
    state = sampler.run_mcmc(p0,Nmc)

## run MCMC with on cpu
#sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
#state = sampler.run_mcmc(p0,Nmc)
    
print("MCMC Total Run Time [s] =", time.time()-t0full)


# # Pickle MCMC Chain
del(sampler.pool)
with open('./output/'+galaxy+'_mcmc.pkl', 'wb') as f:
    pickle.dump(sampler, f, pickle.HIGHEST_PROTOCOL)

# # Unpickle MCMC Chain
with open('./output/'+galaxy+'_mcmc.pkl', 'rb') as f:
    sampler = pickle.load(f)


# # Find best-fit model (max logP) and evaluate

samples = sampler.chain
#print(np.shape(samples))
#print(np.median(samples[:,-1,:],axis=0))

logp=sampler.lnprobability
maxlogp=np.max(logp)
selbest=np.where(logp == maxlogp)
bestwalk=selbest[0][0]
bestsamp=selbest[1][0]
pbest=samples[bestwalk, bestsamp, :]
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
plt.savefig('./plots/'+galaxy+'_corr_model.png')
#plt.show()
        

## Make MCMC Samples Plots

fig, axes = plt.subplots(8, figsize=(20, 50), sharex=True)
labels = ["l", "rc", "tc", "ts", "tfb", "Ng", "voff"]
for i in range(ndim):
    ax = axes[i]
    for j in range(nwalkers):
        ax.plot(samples[j, :, i], "-", alpha=0.5)
#    ax.set_xlim(0, Nmc)
#    ax.set_ylim(-20,100)
    ax.set_ylabel(labels[i], fontsize=20)
#    ax.yaxis.set_label_coords(-0.1, 0.5)
    ax.set_xlabel("step number", fontsize=20)   
ax=axes[7]
for j in range(nwalkers):
    ax.plot(logp[j,:])
    ax.set_ylabel("logP", fontsize=20)
#    ax.set_ylim(100,150)
plt.savefig('./plots/'+galaxy+'_mcmc_samples.png')

## Make MCMC Corner Plot
Nburn=250
goodsamples=samples[:,Nburn:-1,:]
flat_goodsamples=goodsamples.reshape((np.shape(goodsamples)[0]*np.shape(goodsamples)[1],np.shape(goodsamples)[2]))
fig = corner.corner(flat_goodsamples, labels=labels, bins=10, hist_bin_factor=1, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})

plt.savefig('./plots/'+galaxy+'_mcmc_corner.png')
