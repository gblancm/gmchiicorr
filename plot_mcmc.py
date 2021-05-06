import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle5 as pickle
import corner
import emcee


galaxy=sys.argv[1]
print(galaxy)

def log_prob(p):
        return 0

# # Unpickle MCMC Chain
with open('./output/'+galaxy+'_mcmc.pkl', 'rb') as f:
    sampler = pickle.load(f)


## Plot histogram of acceptance fraction
fig, ax = plt.subplots(figsize=(10, 10), sharex=True)
ax.hist(sampler.acceptance_fraction)
ax.set_xlabel('MCMC Acceptance Fraction', fontsize=20)
plt.savefig('./plots/'+galaxy+'_mcmc_acceptance_plot.png')



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


# # Make MCMC Plots

print(np.shape(samples))
nwalkers=np.shape(samples)[0]
ndim=np.shape(samples)[2]

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

plt.savefig('./plots/'+galaxy+'_mcmc_samples_plot.png')


Nburn=500
accept=0
goodsamples=samples[(sampler.acceptance_fraction>=accept),Nburn:-1,:]
flat_goodsamples=goodsamples.reshape((np.shape(goodsamples)[0]*np.shape(goodsamples)[1],np.shape(goodsamples)[2]))
#fig = corner.corner(flat_goodsamples, labels=labels, range=[(100,300), (5,100), (1,500), (1,10), (1,30), (1, 10), (0,30)])
fig = corner.corner(flat_goodsamples, labels=labels, bins=20, hist_bin_factor=1, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})

plt.savefig('./plots/'+galaxy+'_mcmc_corner_plot.png')
