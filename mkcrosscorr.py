## This code measure the 2-point Cross-Correlation function of GMCs and HII Regions in PHANGS galaxies (ALMA+MUSE)

import sys
import numpy as np
import matplotlib.pyplot as plt
from reproject import reproject_interp
from scipy.interpolate import NearestNDInterpolator
from scipy.optimize import curve_fit
import astropy.io.ascii as ascii
import astropy.io.fits as fits
from astropy import wcs




#from scipy.spatial.distance import cdist
#from astropy.coordinates import SkyCoord
#from photutils import SkyCircularAperture
#from photutils import aperture_photometry
#import astropy.units as u
#import time


from gmchiicorr import getparams
from gmchiicorr import radec2xy
from gmchiicorr import rancat
from gmchiicorr import w
from gmchiicorr import egratio
from gmchiicorr import lin

np.random.seed(666)
 
# ## Name of galaxy

galaxy=sys.argv[1]
print("Measuring 2-point Cross-Correlation function for: ", galaxy)


# Get position parameters for galaxy
ra0, dec0, D, PA, inc = getparams(galaxy)
print(galaxy)
print('ra0, dec0 = ', ra0, dec0)
print('D [Mpc]= ', D)
print('PA, inc = ', PA, inc)

# Use NED values for some troublesome galaxies
if galaxy=='NGC1365':
    PA=49.5
    inc=41.7
    print('Warning!!! - Overriding Sample Table Geometric Parameters')
    print('ra0, dec0 = ', ra0, dec0)
    print('D = ', D)
    print('PA, inc = ', PA, inc)

if galaxy=='IC5332':
    PA=27.5
    inc=30.9
    print('Warning!!! - Overriding Sample Table Geometric Parameters')
    print('ra0, dec0 = ', ra0, dec0)
    print('D = ', D)
    print('PA, inc = ', PA, inc)



# Get coordinates of GMCs from ALMA and map onto disk coordinates
gmcdata=fits.open('./catalogs/gmc/'+galaxy.lower()+'_co21_native_props.fits')
ra1=gmcdata[1].data['XCTR_DEG']
dec1=gmcdata[1].data['YCTR_DEG']
x1, y1, r1, theta1, dra1, ddec1 = radec2xy(ra1, dec1, ra0, dec0, D, PA, inc)

# Get coordinates of HII regions from MUSE and map onto disk coordinates
hiidata=ascii.read('./catalogs/hiireg/'+galaxy+'_emission_region_cat.txt')
ra2=hiidata['ra_peak'].data
dec2=hiidata['dec_peak'].data
x2, y2, r2, theta2, dra2, ddec2 = radec2xy(ra2, dec2, ra0, dec0, D, PA, inc)


# ## Create joint mask of regions observed in both CO and Ha
print("Creating ALMA, MUSE, and joint masks.")

# read CO map and create mask (NaN are not observed)
comap=fits.open('./maps/co/'+galaxy.lower()+'_12m+7m+tp_co21_strict_mom0.fits')
cowcs=wcs.WCS(comap[0], naxis=2)
nxco=comap[0].header['NAXIS1']
nyco=comap[0].header['NAXIS2']
comask=np.ones((nyco, nxco), dtype=bool)
comask[np.isnan(comap[0].data)]=False

# read Ha map, reproject to CO map and create mask (NaN are not observed)
hamap=fits.open('./maps/ha/'+galaxy+'_IMAGE_FOV_WFI_NB.fits')
hamap2, hafootprint = reproject_interp(hamap[1], cowcs, shape_out=(nyco, nxco))
hamask=np.ones((nyco, nxco), dtype=bool)
hamask[(hamap2==0)+(np.isnan(hamap2))]=False

# merge both masks into total mask and mask edges
mask=np.logical_and(comask, hamask)
mask[0,:]=False
mask[:,0]=False
mask[-1,:]=False
mask[:,-1]=False

# calculate disk corrdinates of mask pixels
x0mask, y0mask = np.meshgrid(np.arange(nxco), np.arange(nyco))
ramask, decmask = cowcs.all_pix2world(x0mask, y0mask, 0)
xmask, ymask, rmask, thetamask, dramask, ddecmask = radec2xy(ramask, decmask, ra0, dec0, D, PA, inc)

# create mask arrays for GMCs and HII regions by interpolating the mask to the GMC and HII region deprojected coords
aux=np.stack((xmask.flatten(), ymask.flatten()), axis=-1)
f=NearestNDInterpolator(aux, mask.flatten())
sel1=f(x1, y1)
sel2=f(x2, y2)

# Write mask to fits file to be used later
hdu1=fits.PrimaryHDU(mask.astype(int))
hdu1.header=comap[0].header
hdu2=fits.PrimaryHDU(xmask)
hdu3=fits.PrimaryHDU(ymask)
hdul=fits.HDUList([hdu1])
hdul.append(hdu2)
hdul.append(hdu3)
hdul.writeto('./output/'+galaxy+'_mask.fits', overwrite=True)

# Make random catalog 

myNr=50000
xr, yr = rancat(xmask, ymask, mask, f=f, Nr=myNr)

# Calculate HII-GMC number ratio
fhg, efhg = egratio(len(x2[sel2]), len(x1[sel1]))
print("HII regio to GMC ratio (error)= ", fhg, efhg)


## Calculate 2-point cross-correlation function between GMCs and HII regions

r0, w0, ew0 = w(x1[sel1], y1[sel1], x2[sel2], y2[sel2], xr, yr)


## Fit and remove large scale (few kpc) correlation using linear model

rmin=500
rmax=1500

# restricted fitting raneg for some galaxies
if galaxy=='NGC0628':
    rmax=1000


sel=(r0>=rmin)*(r0<=rmax)
popt, pcov = curve_fit(lin, r0[sel], w0[sel], sigma=ew0[sel])
w0small=w0-lin(r0, *popt)



# Write all results to disk
ascii.write([["fhg", "efhg", "Nh", "Ng"], np.array([fhg, efhg, len(x2[sel2]), len(x1[sel1])])], './output/'+galaxy+'_fhg.txt', format='no_header', overwrite=True)
ascii.write([r0, w0, ew0], './output/'+galaxy+'_corr.txt', format='no_header', overwrite=True)
ascii.write([r0[r0<=500], w0small[r0<=500], ew0[r0<=500]], './output/'+galaxy+'_corr_small.txt', format='no_header', overwrite=True)
ascii.write([x1[sel1], y1[sel1]], './output/'+galaxy+'_xy_gmc.txt', format='no_header', overwrite=True)
ascii.write([x2[sel2], y2[sel2]], './output/'+galaxy+'_xy_hii.txt', format='no_header', overwrite=True)
ascii.write([xr, yr], './output/'+galaxy+'_xy_rand.txt', format='no_header', overwrite=True)



# ## Plot GMC - HII region 2-point cross-correlation function


fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(r0, w0, ew0, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0, w0, 'o', color='black', alpha=1.0)
#ax.set_xlim(0, 250)
ax.axhline(y=0, linestyle='--')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'$\omega(r)$ [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
#ax.set_yscale('log')
plt.savefig('./plots/'+galaxy+'_corr.png')
#plt.show()


fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(r0, w0, ew0, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0, w0, 'o', color='black', alpha=1.0)
ax.set_xlim(0, 500)
ax.axhline(y=0, linestyle='--')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'$\omega(r)$ [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
#ax.set_yscale('log')
plt.savefig('./plots/'+galaxy+'_corr_zoom.png')
#plt.show()

fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(r0, w0, ew0, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0, w0, 'o', color='black', alpha=1.0)
#ax.set_xlim(0, 250)
ax.set_ylim(1e-3, 1e1)
ax.axhline(y=0, linestyle='--')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'log($\omega(r)$) [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
ax.set_yscale('log')
plt.savefig('./plots/'+galaxy+'_corr_log.png')
#plt.show()

#Plot large scale fit and small scale cross-correlation function

fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(r0, w0, ew0, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0, w0, 'o', color='black', alpha=1.0)
ax.plot(r0[sel], w0[sel], 'o', color='red', alpha=0.5)
ax.plot(r0, lin(r0, *popt), color='red')
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
ax.errorbar(r0, w0small, ew0, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0, w0small, 'o', color='black', alpha=1.0)
ax.plot(r0[sel], w0small[sel], 'o', color='red', alpha=0.5)
ax.set_xlim(0, 500)
#ax.set_ylim(1e-3, 2e0)
ax.axhline(y=0, linestyle='--')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'$\omega(r)-\omega_{lin}$ [pc]', fontsize=20)
plt.title(galaxy+" ; fhg="+"{:.2f}".format(fhg)+" ("+"{:.2f}".format(efhg)+")", fontsize=30)
ax.tick_params(labelsize=20)
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.savefig('./plots/'+galaxy+'_corr_small.png')
#plt.show()




# ## Plot GMCs with HII regions and Random Catalog


fig, ax = plt.subplots(figsize=(12, 12))
ax.plot(x1[sel1], y1[sel1], 'ob', alpha=0.3, label='GMCs')
ax.plot(x2[sel2], y2[sel2], 'or', alpha=0.3, label='HII Regions')
ax.plot(x1, y1, '.b', alpha=0.3)
ax.plot(x2, y2, '.r', alpha=0.3)
#ax.plot(xr, yr, '.g', alpha=0.1)
ax.set_xlim(np.min(xmask[mask]), np.max(xmask[mask]))
ax.set_ylim(np.min(ymask[mask]), np.max(ymask[mask]))
ax.contour(xmask, ymask, mask)
ax.set_aspect('equal', 'datalim')
ax.axis('equal')
plt.xlabel('X [pc]', fontsize=20)
plt.ylabel('Y [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)
plt.savefig('./plots/'+galaxy+'_xy.png')
#plt.show()

fig, ax = plt.subplots(figsize=(12, 12))
ax.plot(x1[sel1], y1[sel1], 'ob', alpha=0.3, label='GMCs')
#ax.plot(x2[sel2], y2[sel2], 'or', alpha=0.3)
ax.plot(x1, y1, '.b', alpha=0.3)
#ax.plot(x2, y2, '.r', alpha=0.3)
ax.plot(xr, yr, '.g', alpha=0.1, label='Random Points')
ax.set_xlim(np.min(xmask[mask]), np.max(xmask[mask]))
ax.set_ylim(np.min(ymask[mask]), np.max(ymask[mask]))
ax.contour(xmask, ymask, mask)
ax.set_aspect('equal', 'datalim')
ax.axis('equal')
plt.xlabel('X [pc]', fontsize=20)
plt.ylabel('Y [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)
plt.savefig('./plots/'+galaxy+'_xy_rand.png')
#plt.show()

fig, ax = plt.subplots(figsize=(12, 12))
ax.plot(x1[sel1], y1[sel1], 'ob', alpha=0.3, label='GMCs')
ax.plot(x2[sel2], y2[sel2], 'or', alpha=0.3, label='HII Regions')
ax.plot(x1, y1, '.b', alpha=0.3)
ax.plot(x2, y2, '.r', alpha=0.3)
#ax.plot(xr, yr, '.g', alpha=0.1)
ax.set_xlim(-1000, 1000)
ax.set_ylim(-1000, 1000)
ax.contour(xmask, ymask, mask)
plt.xlabel('X [pc]', fontsize=20)
plt.ylabel('Y [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
ax.legend(fontsize=20, loc=1)
plt.savefig('./plots/'+galaxy+'_xy_zoom.png')
#plt.show()

fig, ax = plt.subplots(figsize=(12, 12))
ax.plot(x1[sel1], y1[sel1], 'ob', alpha=0.3, label='GMCs')
#ax.plot(x2[sel2], y2[sel2], 'or', alpha=0.3)
ax.plot(x1, y1, '.b', alpha=0.3)
#ax.plot(x2, y2, '.r', alpha=0.3)
ax.plot(xr, yr, '.g', alpha=0.1, label='Random Points')
ax.set_xlim(-1000, 1000)
ax.set_ylim(-1000, 1000)
ax.contour(xmask, ymask, mask)
plt.xlabel('X [pc]', fontsize=20)
plt.ylabel('Y [pc]', fontsize=20)
plt.title(galaxy, fontsize=30)
ax.tick_params(labelsize=20)
ax.legend(fontsize=20)
plt.savefig('./plots/'+galaxy+'_xy_rand_zoom.png')
#plt.show()




