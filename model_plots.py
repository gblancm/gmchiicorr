## This code makes plots of the cross-correlation function model changing one parameter at a time around a fiducial model

import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.ascii as ascii

from gmchiicorr import eval_w


# test bins feature in eval_w
#bins=np.array([15, 45, 75])
bins=None

# Define fiducial model
rmax=300
pbest=[200, 75, 10, 5, 2, 10, 15]
r00, w00, ew00, fhg00 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins)


# Plots for "voff"
parr=np.array([0, 5, 10, 15, 20])
fig, ax = plt.subplots(figsize=(12, 8))
infiles=glob.glob('./output/*corr_small.txt')
for i in range(len(infiles)):
    tab=ascii.read(infiles[i])
    auxr=tab['col1'].data
    auxerr=tab['col2'].data
    ax.plot(auxr, auxerr, ':', color='grey', alpha=0.3)
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Fiducial ; fhg='+"{:.2f}".format(fhg00))    
for i in range(len(parr)):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=parr[i], rmax=rmax, bins=bins)
    ax.plot(r0, w0, '-o', alpha=1.0, label='voff='+str(parr[i])+' ; fhg='+"{:.2f}".format(fhg0))    
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/model_plots_voff.png')
#plt.show()

# Plots for "Ng"
parr=np.array([1, 2, 5, 10, 20])
fig, ax = plt.subplots(figsize=(12, 8))
for i in range(len(infiles)):
    tab=ascii.read(infiles[i])
    auxr=tab['col1'].data
    auxerr=tab['col2'].data
    ax.plot(auxr, auxerr, ':', color='grey', alpha=0.3)
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Fiducial ; fhg='+"{:.2f}".format(fhg00))    
for i in range(len(parr)):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=parr[i], voff0=pbest[6], rmax=rmax, bins=bins)
    ax.plot(r0, w0, '-o', alpha=1.0, label='Ng='+str(parr[i])+' ; fhg='+"{:.2f}".format(fhg0))    
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/model_plots_Ng.png')
#plt.show()


# Plots for "tfb"
parr=np.array([1, 2, 3, 4, 5])
fig, ax = plt.subplots(figsize=(12, 8))
for i in range(len(infiles)):
    tab=ascii.read(infiles[i])
    auxr=tab['col1'].data
    auxerr=tab['col2'].data
    ax.plot(auxr, auxerr, ':', color='grey', alpha=0.3)
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Fiducial ; fhg='+"{:.2f}".format(fhg00))    
for i in range(len(parr)):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=parr[i], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins)
    ax.plot(r0, w0, '-o', alpha=1.0, label='tfb='+str(parr[i])+' ; fhg='+"{:.2f}".format(fhg0))    
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/model_plots_tfb.png')
#plt.show()


# Plots for "ts"
parr=np.array([1, 2.5, 5.0, 7.5, 10])
fig, ax = plt.subplots(figsize=(12, 8))
for i in range(len(infiles)):
    tab=ascii.read(infiles[i])
    auxr=tab['col1'].data
    auxerr=tab['col2'].data
    ax.plot(auxr, auxerr, ':', color='grey', alpha=0.3)
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Fiducial ; fhg='+"{:.2f}".format(fhg00))    
for i in range(len(parr)):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=pbest[2], ts0=parr[i], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins)
    ax.plot(r0, w0, '-o', alpha=1.0, label='ts='+str(parr[i])+' ; fhg='+"{:.2f}".format(fhg0))    
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/model_plots_ts.png')
#plt.show()


# Plots for "tc"
parr=np.array([1, 3, 5, 10, 30])
fig, ax = plt.subplots(figsize=(12, 8))
for i in range(len(infiles)):
    tab=ascii.read(infiles[i])
    auxr=tab['col1'].data
    auxerr=tab['col2'].data
    ax.plot(auxr, auxerr, ':', color='grey', alpha=0.3)
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Fiducial ; fhg='+"{:.2f}".format(fhg00))    
for i in range(len(parr)):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=pbest[1], tc0=parr[i], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins)
    ax.plot(r0, w0, '-o', alpha=1.0, label='tc='+str(parr[i])+' ; fhg='+"{:.2f}".format(fhg0))    
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/model_plots_tc.png')
#plt.show()


# Plots for "rc"
parr=np.array([25, 50, 75, 100, 125])
fig, ax = plt.subplots(figsize=(12, 8))
for i in range(len(infiles)):
    tab=ascii.read(infiles[i])
    auxr=tab['col1'].data
    auxerr=tab['col2'].data
    ax.plot(auxr, auxerr, ':', color='grey', alpha=0.3)
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Fiducial ; fhg='+"{:.2f}".format(fhg00))    
for i in range(len(parr)):
    r0, w0, ew0, fhg0 = eval_w(l0=pbest[0], rc0=parr[i], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins)
    ax.plot(r0, w0, '-o', alpha=1.0, label='rc='+str(parr[i])+' ; fhg='+"{:.2f}".format(fhg0))
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/model_plots_rc.png')
#plt.show()


# Plots for "l"
parr=np.array([50, 100, 150, 200, 250])
fig, ax = plt.subplots(figsize=(12, 8))
for i in range(len(infiles)):
    tab=ascii.read(infiles[i])
    auxr=tab['col1'].data
    auxerr=tab['col2'].data
    ax.plot(auxr, auxerr, ':', color='grey', alpha=0.3)
ax.plot(r00, w00, '-o', alpha=1.0, color='black', label='Fiducial ; fhg='+"{:.2f}".format(fhg00))    
for i in range(len(parr)):
    r0, w0, ew0, fhg0 = eval_w(l0=parr[i], rc0=pbest[1], tc0=pbest[2], ts0=pbest[3], tfb0=pbest[4], Ng0=pbest[5], voff0=pbest[6], rmax=rmax, bins=bins)
    ax.plot(r0, w0, '-o', alpha=1.0, label='l='+str(parr[i])+' ; fhg='+"{:.2f}".format(fhg0))  
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle='--')
ax.set_xlabel('r [pc]', fontsize=20)
ax.set_ylabel(r'$\omega(r)$ [pc]', fontsize=20)
ax.legend(fontsize=20)
ax.tick_params(labelsize=20)
ax.set_title('Fiducial: l='+str(pbest[0])+'; rc='+str(pbest[1])+'; tc='+str(pbest[2])+'; ts='+str(pbest[3])+'; tfb='+str(pbest[4])+'; Ng='+str(pbest[5])+'; voff='+str(pbest[6]), fontsize=15)
plt.savefig('./plots/model_plots_l.png')
#plt.show()
