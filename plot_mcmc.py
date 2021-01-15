import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle5 as pickle
import corner
import emcee



galaxy=sys.argv[1]
print(galaxy)


def log_prob(p):
    lprior=log_prior(p)
    if np.isfinite(lprior):
        r0, w0, ew0 = eval_w(l0=p[0], rc0=p[1], tc0=p[2], ts0=p[3], tfb0=p[4], Ng0=p[5], voff0=p[6])
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


Nburn=250
goodsamples=samples[:,Nburn:-1,:]
flat_goodsamples=goodsamples.reshape((np.shape(goodsamples)[0]*np.shape(goodsamples)[1],np.shape(goodsamples)[2]))
#fig = corner.corner(flat_goodsamples, labels=labels, range=[(100,300), (5,100), (1,500), (1,10), (1,30), (1, 10), (0,30)])
fig = corner.corner(flat_goodsamples, labels=labels)

plt.savefig('./plots/'+galaxy+'_mcmc_corner_plot.png')
