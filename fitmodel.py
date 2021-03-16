import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
import astropy.io.ascii as ascii
import time
from scipy.optimize import curve_fit
import emcee
import pickle
import os
from multiprocessing import Pool
from multiprocessing import cpu_count
os.environ["OMP_NUM_THREADS"] = "1"
import corner
from fastdist import fastdist
from fast_histogram import histogram1d


from gmchiicorr import w
from gmchiicorr import drawgmc_l
from gmchiicorr import drawhii
from gmchiicorr import lin
from gmchiicorr import eval_w


#np.random.seed(666)

# ## Name of galaxy

#galaxy='IC5332'
#galaxy='NGC0628'
#galaxy='NGC1087'
#galaxy='NGC1365' 
#galaxy='NGC1512'
#galaxy='NGC1672'
#galaxy='NGC2835'
#galaxy='NGC3351'
#galaxy='NGC3627'
#galaxy='NGC5068'
galaxy=sys.argv[1]
print(galaxy)


# ## Read GMC and Random catalogs coordinates

xygmc=ascii.read('./output/'+galaxy+'_xy_gmc.txt')
x1=xygmc['col0'].data
y1=xygmc['col1'].data

xyrand=ascii.read('./output/'+galaxy+'_xy_rand.txt')
xr=xyrand['col0'].data
yr=xyrand['col1'].data


# ## Read Observed Correlation Function

obscorr=ascii.read('./output/'+galaxy+'_corr.txt')
r0obs=obscorr['col0'].data
w0obs=obscorr['col1'].data
ew0obs=obscorr['col2'].data


# ## Fit and remove large scale (few kpc) correlation using linear model

rmin=500
rmax=3000
sel=(r0obs>=rmin)*(r0obs<=rmax)
popt, pcov = curve_fit(lin, r0obs[sel], w0obs[sel], sigma=ew0obs[sel])
#print(popt)
w0small=w0obs-lin(r0obs, *popt)


#plot large scale fit

fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(r0obs, w0obs, ew0obs, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0obs, w0obs, 'o', color='black', alpha=1.0)
ax.plot(r0obs[sel], w0obs[sel], 'o', color='red', alpha=0.5)
ax.plot(r0obs, lin(r0obs, *popt), color='red')
#ax.set_ylim(1e-3, 2e0)
ax.axhline(y=0, linestyle='--')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'$\omega(r)$ [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.savefig('./plots/'+galaxy+'_corr_linfit.png')
#plt.show()

fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(r0obs, w0small, ew0obs, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0obs, w0small, 'o', color='black', alpha=1.0)
ax.plot(r0obs[sel], w0small[sel], 'o', color='red', alpha=0.5)
ax.set_xlim(0, 1e3)
#ax.set_ylim(1e-3, 2e0)
ax.axhline(y=0, linestyle='--')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'$\omega(r)-\omega_{lin}$ [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.savefig('./plots/'+galaxy+'_corr_small.png')
#plt.show()



# # Define Priors and Likelihood Functions

# Trim observed corr function to <=1kpc
selr=(r0obs<=1000)

lrange=np.array([50,500])
rcrange=np.array([5,100])
tcrange=np.array([1,500])
#tsrange=np.array([1,20])
tsrange=np.array([3,5])  # from Chevance et al.
#tfbrange=np.array([1,20])
tfbrange=np.array([0.1,5]) # for Chevance et al. ts prior
Ngrange=np.array([1,10])
voffrange=np.array([0,30])

def log_prior(p):
    l1, rc1, tc1, ts1, tfb1, Ng1, voff1 = p
    if lrange[0]<l1<lrange[1] and rcrange[0]<rc1<rcrange[1] and tcrange[0]<tc1<tcrange[1] and tsrange[0]<ts1<tsrange[1] and tfbrange[0]<tfb1<tfbrange[1] and Ngrange[0]<Ng1<Ngrange[1] and voffrange[0]<voff1<voffrange[1]:
        return 0.0
    else:
        return -np.inf



def log_prob(p):
    lprior=log_prior(p)
    if np.isfinite(lprior):
        r0, w0, ew0, fhg0 = eval_w(l0=p[0], rc0=p[1], tc0=p[2], ts0=p[3], tfb0=p[4], Ng0=p[5], voff0=p[6])
        res=w0-w0small[selr]
        sig=ew0obs[selr]
        prob=1/(2*np.pi*sig**2)*np.exp(-0.5*(res/sig)**2)
        logp=lprior+np.sum(np.log(prob))
        if not np.isfinite(logp):
#        print("fail2")
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


with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
    state = sampler.run_mcmc(p0,Nmc)

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
print(np.shape(samples))
print(np.median(samples[:,-1,:],axis=0))

logp=sampler.lnprobability
maxlogp=np.max(logp)
selbest=np.where(logp == maxlogp)
bestwalk=selbest[0][0]
bestsamp=selbest[1][0]
pbest=samples[bestwalk, bestsamp, :]
print("Best-fit Parameters:", pbest)

r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6])

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(r0, w0, 'o', color='green', alpha=1.0)    
ax.errorbar(r0obs, w0small, ew0obs, fmt="o", color='black', capsize=5, alpha=0.5)
ax.plot(r0obs, w0small, '-o', color='black', alpha=0.5)
ax.set_xlim(0, 1000)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.set_title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
plt.savefig('./plots/'+galaxy+'_corr_model.png')
#plt.show()
        
#fig, ax = plt.subplots(figsize=(12, 12))
##ax.plot(xgmc, ygmc, '.', color='blue', alpha=0.4)
#ax.plot(xgmc[fgmc], ygmc[fgmc], 'o', color='blue', label='GMC')
##ax.plot(xhii, yhii, '.', color='red', alpha=0.4)
#ax.plot(xhii[fhii], yhii[fhii], 'o', color='red', label='HII')
##ax.set_xlim(-1000, 1000)
##ax.set_ylim(-1000, 1000)
#ax.set_xlabel('X [pc]', fontsize=20)
#ax.set_ylabel('Y [pc]', fontsize=20)
#ax.set_title(galaxy, fontsize=30)
#ax.tick_params(labelsize=20)
#ax.legend(fontsize=20)
#plt.savefig('./plots/'+galaxy+'_xy_model.png')
##plt.show()


# # Make MCMC Plots


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


Nburn=0
goodsamples=samples[:,Nburn:-1,:]
flat_goodsamples=goodsamples.reshape((np.shape(goodsamples)[0]*np.shape(goodsamples)[1],np.shape(goodsamples)[2]))
fig = corner.corner(flat_goodsamples, labels=labels, range=[(100,300), (5,100), (1,500), (1,10), (1,30), (1, 10), (0,30)])

plt.savefig('./plots/'+galaxy+'_mcmc_corner.png')
