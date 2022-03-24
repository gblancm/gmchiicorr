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


# Find initial guess via simple minimization
print("====================================")
print("Finding Starting Point via curve_fit")
print("====================================")

## All free parameters
#def func1(r, p0, p1, p2, p3, p4, p5, p6):
#    bins=r[0:-1]
#    r0, w0, ew0, fhg0 = eval_w(l0=p0, rc0=p1, tc0=p2, ts0=p3, tfb0=p4, Ng0=p5, voff0=p6, bins=bins, Nsamples=150)  #Nsamples=150 yields rms smaller than measurement errors
#    return np.concatenate([w0,np.array([fhg0])])    
#p0=np.array([200, 100, 10, 5, 2, 5, 5])
#auxr=np.concatenate([r0obs,np.array([-1])])
#auxw=np.concatenate([w0obs,np.array([fhgobs])])
#auxew=np.concatenate([ew0obs,np.array([efhgobs])])
#pstart, pcov = curve_fit(func1, auxr, auxw, p0=p0, sigma=auxew, method='lm', epsfcn=0.01)

### Fixing ts=5
#def func1(r, p0, p1, p2, p4, p5, p6):
#    bins=r[0:-1]
#    r0, w0, ew0, fhg0 = eval_w(l0=p0, rc0=p1, tc0=p2, ts0=5, tfb0=p4, Ng0=p5, voff0=p6, bins=bins, Nsamples=150)  #Nsamples=150 yields rms smaller than measurement errors
#    return np.concatenate([w0,np.array([fhg0])])    
#p0=np.array([200, 100, 10, 2, 5, 5])
#auxr=np.concatenate([r0obs,np.array([-1])])
#auxw=np.concatenate([w0obs,np.array([fhgobs])])
#auxew=np.concatenate([ew0obs,np.array([efhgobs])])
#pstart, pcov = curve_fit(func1, auxr, auxw, p0=p0, sigma=auxew, method='lm', epsfcn=0.01)
#pstart=np.array([pstart[0], pstart[1], pstart[2], 5.0, pstart[3], pstart[4], pstart[5]])

## Fixing ts=5 and Ng=10
def func1(r, p0, p1, p2, p4, p6):
    bins=r[0:-1]
    r0, w0, ew0, fhg0 = eval_w(l0=p0, rc0=p1, tc0=p2, ts0=5, tfb0=p4, Ng0=10, voff0=p6, bins=bins, Nsamples=150)  #Nsamples=150 yields rms smaller than measurement errors
    return np.concatenate([w0,np.array([fhg0])])    
p0=np.array([200, 100, 10, 2, 5])
auxr=np.concatenate([r0obs,np.array([-1])])
auxw=np.concatenate([w0obs,np.array([fhgobs])])
auxew=np.concatenate([ew0obs,np.array([efhgobs])])
pstart, pcov = curve_fit(func1, auxr, auxw, p0=p0, sigma=auxew, method='lm', epsfcn=0.01)
pstart=np.array([pstart[0], pstart[1], pstart[2], 5.0, pstart[3], 10.0, pstart[4]])




## Plot pstart fit
pbest=pstart
print("Best-fit Parameters:", pbest)
r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], bins=r0obs)
fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(r0, w0, '-o', color='green', alpha=1.0, label="fhg="+"{:.2f}".format(fhg0))    
ax.errorbar(r0obs, w0obs, ew0obs, fmt="o", color='black', capsize=5, alpha=0.5)
ax.plot(r0obs, w0obs, 'o', color='black', alpha=0.5)
#ax.set_xlim(0, 00)
ax.axhline(y=0, linestyle=':', color='grey')

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
plt.savefig('./plots/'+galaxy+'_pstart.png')


# Define prior parameter space
lrange=np.array([50,500])
rcrange=np.array([5,300])
tcrange=np.array([1,50])
tsrange=np.array([0,15])
tfbrange=np.array([0,15])
Ngrange=np.array([1,30])
voffrange=np.array([0,50])

# Define uniform prior distribution
def log_prior(p):
    l1, rc1, tc1, ts1, tfb1, Ng1, voff1 = p
    if lrange[0]<l1<lrange[1] and rcrange[0]<rc1<rcrange[1] and tcrange[0]<tc1<tcrange[1] and tsrange[0]<ts1<tsrange[1] and tfbrange[0]<tfb1<tfbrange[1] and Ngrange[0]<Ng1<Ngrange[1] and voffrange[0]<voff1<voffrange[1]:
        return 0.0
    else:
        return -np.inf


# Define normal prior for ts1
def log_tsprior(p):
    l1, rc1, tc1, ts1, tfb1, Ng1, voff1 = p
    mu = 5
    sigma = 2
    return np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(ts1-mu)**2/sigma**2

# Define normal prior for Ng
def log_ngprior(p):
    l1, rc1, tc1, ts1, tfb1, Ng1, voff1 = p
    mu = 10
    sigma = 5
    return np.log(1.0/(np.sqrt(2*np.pi)*sigma))-0.5*(Ng1-mu)**2/sigma**2


# Define likelihood*prior distriibution
def log_prob(p):
    #lprior=log_prior(p)  # without extra prior in Ng
    lprior=log_prior(p)+log_tsprior(p)+log_ngprior(p)    # for Gaussian Prior on ts and Ng
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


## Set up MCMC
ndim=7
#nwalkers=256
#Nmc=3000
nwalkers=32
Nmc=200

p0 = np.zeros((nwalkers, ndim))
## Initialize walkers uniformly across prior paramter space
#p0[:,0]=np.random.uniform(lrange[0], lrange[1], nwalkers)
#p0[:,1]=np.random.uniform(rcrange[0], rcrange[1], nwalkers)
#p0[:,2]=np.random.uniform(tcrange[0], tcrange[1], nwalkers)
#p0[:,3]=np.random.uniform(tsrange[0], tsrange[1], nwalkers)
#p0[:,4]=np.random.uniform(tfbrange[0], tfbrange[1], nwalkers)
#p0[:,5]=np.random.uniform(Ngrange[0], Ngrange[1], nwalkers)
#p0[:,6]=np.random.uniform(voffrange[0], voffrange[1], nwalkers)

## Initialize walkers in a cloud around pstart (from curve_fit)
factor=0.25
p0[:,0]=pstart[0]+np.random.normal(0, factor*pstart[0], nwalkers)
p0[:,1]=pstart[1]+np.random.normal(0, factor*pstart[1], nwalkers)
p0[:,2]=pstart[2]+np.random.normal(0, factor*pstart[2], nwalkers)
p0[:,3]=pstart[3]+np.random.normal(0, factor*pstart[3], nwalkers)
p0[:,4]=pstart[4]+np.random.normal(0, factor*pstart[4], nwalkers)
p0[:,5]=pstart[5]+np.random.normal(0, factor*pstart[5], nwalkers)
p0[:,6]=pstart[6]+np.random.normal(0, factor*pstart[6], nwalkers)



# # Run MCMC Chain
print("====================================")
print("Starting MCMC")
print("====================================")

t0full=time.time()

# run MCMC with multiple cpus
with Pool() as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool, moves=[(emcee.moves.StretchMove(a=1.2), 1.0), (emcee.moves.DEMove(), 0.0),(emcee.moves.DESnookerMove(), 0.0),])
    state = sampler.run_mcmc(p0,Nmc)

### run MCMC with on cpu
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



## Plot histogram of acceptance fraction
fig, ax = plt.subplots(figsize=(10, 10), sharex=True)
ax.hist(sampler.acceptance_fraction)
ax.set_xlabel('MCMC Acceptance Fraction', fontsize=20)
plt.savefig('./plots/'+galaxy+'_mcmc_acceptance.png')


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
ax.axhline(y=0, linestyle=':', color='grey')
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
Nburn=0
goodsamples=samples[:,Nburn:-1,:]
flat_goodsamples=goodsamples.reshape((np.shape(goodsamples)[0]*np.shape(goodsamples)[1],np.shape(goodsamples)[2]))
fig = corner.corner(flat_goodsamples, labels=labels, bins=10, hist_bin_factor=1, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})

plt.savefig('./plots/'+galaxy+'_mcmc_corner.png')
