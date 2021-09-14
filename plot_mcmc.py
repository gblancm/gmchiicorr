import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle5 as pickle
import corner
import astropy.io.ascii as ascii
from gmchiicorr import eval_w



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


# # Find best-fit model (max logP)

samples = sampler.chain
nwalkers=np.shape(samples)[0]
ndim=np.shape(samples)[2]

print(np.shape(samples))
print(np.median(samples[:,-1,:],axis=0))

logp=sampler.lnprobability
maxlogp=np.max(logp)
selbest=np.where(logp == maxlogp)
bestwalk=selbest[0][0]
bestsamp=selbest[1][0]
pbest=samples[bestwalk, bestsamp, :]
print("Best-fit Parameters:", pbest)


## Find best-fitr model (median of marginalized PDFs)


# ## Read Observed Correlation Function
obscorr=ascii.read('./output/'+galaxy+'_corr_small.txt')
r0obs=obscorr['col1'].data
w0obs=obscorr['col2'].data
ew0obs=obscorr['col3'].data
## Read fraction of HII regions to GMCs
fhgtab=ascii.read('./output/'+galaxy+'_fhg.txt')
fhgobs=fhgtab['col2'].data[0]   
efhgobs=fhgtab['col2'].data[1]

# Get medians of PDF
Nburn=1000 # burn period
pbest=np.zeros(ndim)
for i in range(ndim):
    pbest[i]=np.array(np.median(samples[:,Nburn:-1,i]))
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
plt.savefig('./plots/'+galaxy+'_corr_best_plot.png')
#plt.show()
        



# # Make MCMC Plots

print(np.shape(samples))

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


accept=0.0 # minimum acceptance fraction
goodsamples=samples[(sampler.acceptance_fraction>=accept),Nburn:-1,:]
flat_goodsamples=goodsamples.reshape((np.shape(goodsamples)[0]*np.shape(goodsamples)[1],np.shape(goodsamples)[2]))
#fig = corner.corner(flat_goodsamples, labels=labels, range=[(100,300), (5,100), (1,500), (1,10), (1,30), (1, 10), (0,30)])
fig = corner.corner(flat_goodsamples, labels=labels, bins=20, hist_bin_factor=1, quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})

plt.savefig('./plots/'+galaxy+'_mcmc_corner_plot.png')





