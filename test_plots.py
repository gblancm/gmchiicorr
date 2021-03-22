## This code makes plots of the cross-correlation function model changing one parameter at a time around a fiducial model

import numpy as np
import matplotlib.pyplot as plt

from gmchiicorr import eval_w


# test bins feature in eval_w
#bins=np.array([15, 45, 75])
bins=None

# Define fiducial model
rmax=500
pbest=[200, 50, 20, 10, 2, 5, 10]
#pbest=[247, 75, 10, 3.8, 0.6, 10, 16.2]
r00, w00, ew00, fhg00 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins, Nsamples=500)


# Plots for Nsamples=???
Nsamples=500

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Ns='+"{:.0f}".format(Nsamples)+' ; fhg='+"{:.2f}".format(fhg00))    
for i in range(10):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins, Nsamples=Nsamples)
    ax.plot(r0, w0, '-o', alpha=1.0)    
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/test_plots_Ns'+"{:.0f}".format(Nsamples)+'.png')
#plt.show()

fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(r00, (w00-w00)/w00, '-o', alpha=1.0, color='black', label='Ns='+"{:.0f}".format(Nsamples)+' ; fhg='+"{:.2f}".format(fhg00))    
for i in range(10):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins, Nsamples=Nsamples)
    ax.plot(r0, (w0-w00)/w00, '-o', alpha=1.0)    
ax.set_xlim(0, rmax)
ax.set_ylim(-0.5, 0.5)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\Delta \omega(r) / \omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/test_plots_Ns'+"{:.0f}".format(Nsamples)+'_res.png')
#plt.show()


