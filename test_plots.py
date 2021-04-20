## This code makes plots of the cross-correlation function model changing one parameter at a time around a fiducial model

import numpy as np
import matplotlib.pyplot as plt
from gmchiicorr import eval_w
import glob
import astropy.io.ascii as ascii



# test bins feature in eval_w
#bins=np.array([15, 45, 75])
bins=None

# Define fiducial model
rmax=300
pbest=[200, 50, 20, 10, 2, 5, 10]
#pbest=[247, 75, 10, 3.8, 0.6, 10, 16.2]
r00, w00, ew00, fhg00 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins, Nsamples=1000)


# Plots for Nsamples=???
Nsamples=150




fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(r00, (w00-w00), '-o', alpha=1.0, color='black', label='Ns='+"{:.0f}".format(Nsamples)+' ; fhg='+"{:.2f}".format(fhg00))    
for i in range(30):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins, Nsamples=Nsamples)
    ax.plot(r0, (w0-w00), '-o', alpha=1.0)    
# Plot errors in *smal_corr.txt files
infiles=glob.glob('./output/*corr_small.txt')
for i in range(len(infiles)):
    tab=ascii.read(infiles[i])
    auxr=tab['col1'].data
    auxerr=tab['col3'].data
    ax.plot(auxr, auxerr, ':', color='grey', alpha=0.3)
ax.set_xlim(0, rmax)
ax.set_ylim(-0.5, 0.5)
ax.axhline(0, linestyle='--')
ax.axhline(0.1, linestyle=':')
ax.axhline(-0.1, linestyle=':')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\Delta \omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/test_plots_Ns'+"{:.0f}".format(Nsamples)+'_res.png')
#plt.show()



fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Ns=1000 ; fhg='+"{:.2f}".format(fhg00))    
for i in range(5):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins, Nsamples=Nsamples)
    ax.plot(r0, w0, '-o', alpha=1.0, label='Ns='+"{:.0f}".format(Nsamples)+' ; fhg='+"{:.2f}".format(fhg0))    
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/test_plots_Ns'+"{:.0f}".format(Nsamples)+'.png')
#plt.show()

