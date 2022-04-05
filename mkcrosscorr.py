## This code measure the 2-point Cross-Correlation function of GMCs and HII Regions in PHANGS galaxies (ALMA+MUSE)

import sys
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from reproject import reproject_interp
from scipy.interpolate import NearestNDInterpolator
from scipy.optimize import curve_fit
import astropy.io.ascii as ascii
import astropy.io.fits as fits
from astropy import wcs
import glob
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D

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
#gmcdata=fits.open('./catalogs/gmc/'+galaxy.lower()+'_co21_native_props.fits')
gmcdata=fits.open('./catalogs/v4p0_gmccats/homogen_120pc_148mK/'+galaxy.lower()+'_12m+7m+tp_co21_120pc_148mK_props.fits')
ra1=gmcdata[1].data['XCTR_DEG']
dec1=gmcdata[1].data['YCTR_DEG']
mgmc=gmcdata[1].data['MLUM_MSUN']
x1, y1, r1, theta1, dra1, ddec1 = radec2xy(ra1, dec1, ra0, dec0, D, PA, inc)

# Get coordinates of HII regions from MUSE and map onto disk coordinates
#hiidata=ascii.read('./catalogs/hiireg/'+galaxy+'_emission_region_cat.txt')
#ra2=hiidata['ra_peak'].data
#dec2=hiidata['dec_peak'].data
hiitab=fits.open('./catalogs/Nebulae_catalogue_v2.fits')       
hiidata=hiitab[1].data
hiidatagal=hiidata[(hiidata['gal_name']==galaxy)*(hiidata['HII_class']==1)] # only HII regions in the galaxy of interest 
ra2=hiidatagal['cen_ra']
dec2=hiidatagal['cen_dec']
lha=hiidatagal['Lum_HA6562_CORR']
x2, y2, r2, theta2, dra2, ddec2 = radec2xy(ra2, dec2, ra0, dec0, D, PA, inc)


# ## Create joint mask of regions observed in both CO and Ha
print("Creating ALMA, MUSE, and joint masks.")

# read CO map and create mask (NaN are not observed)
comap=fits.open('./maps/co/'+galaxy.lower()+'_12m+7m+tp_co21_broad_mom0.fits')
cowcs=wcs.WCS(comap[0], naxis=2)
nxco=comap[0].header['NAXIS1']
nyco=comap[0].header['NAXIS2']
comask=np.ones((nyco, nxco), dtype=bool)
comask2=np.ones((nyco, nxco), dtype=bool)
comask[(np.isnan(comap[0].data))+(comap[0].data==0)]=False
comask2[(np.isnan(comap[0].data))]=False

##comask=comask2  # TESTING using footprint

# read Ha map, reproject to CO map and create mask (NaN are not observed)
hamap=fits.open(glob.glob('./maps/ha/'+galaxy+'_IMAGE_FOV_WFI_NB_WCS_Pall_mad_copt_*.fits')[0])
hamap2, hafootprint = reproject_interp(hamap[1], cowcs, shape_out=(nyco, nxco))
hamask=np.ones((nyco, nxco), dtype=bool)
hamask[(hamap2==0)+(np.isnan(hamap2))]=False

# merge both masks into total mask and mask edges
mask=np.logical_and(comask, hamask)
mask[0,:]=False
mask[:,0]=False
mask[-1,:]=False
mask[:,-1]=False

comask2[0,:]=False
comask2[:,0]=False
comask2[-1,:]=False
comask2[:,-1]=False

hamask[0,:]=False
hamask[:,0]=False
hamask[-1,:]=False
hamask[:,-1]=False

# calculate disk corrdinates of mask pixels
x0mask, y0mask = np.meshgrid(np.arange(nxco), np.arange(nyco))
ramask, decmask = cowcs.all_pix2world(x0mask, y0mask, 0)
xmask, ymask, rmask, thetamask, dramask, ddecmask = radec2xy(ramask, decmask, ra0, dec0, D, PA, inc)

# create mask arrays for GMCs and HII regions by interpolating the mask to the GMC and HII region deprojected coords
aux=np.stack((xmask.flatten(), ymask.flatten()), axis=-1)
f=NearestNDInterpolator(aux, mask.flatten())
sel1=f(x1, y1)
sel2=f(x2, y2)

# Print median GMC mass and median Ha luminosity
print("Median GMC Mass = ", np.median(mgmc[sel1]))
print("Median Ha Luminosity = ", np.median(lha[sel2]))


# Write mask to fits file to be used later
hdu1=fits.PrimaryHDU(mask.astype(int))
hdu1.header=comap[0].header
hdu2=fits.PrimaryHDU(xmask)
hdu3=fits.PrimaryHDU(ymask)
hdu4=fits.PrimaryHDU(comask2.astype(int))
hdu5=fits.PrimaryHDU(hamask.astype(int))
hdul=fits.HDUList([hdu1])
hdul.append(hdu2)
hdul.append(hdu3)
hdul.append(hdu4)
hdul.append(hdu5)
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

rmin=200
#rmax=500
rmax=1000

## restricted fitting raneg for some galaxies
#if galaxy=='NGC3351':
#    rmax=1200
#if galaxy=='NGC3627':
#    rmax=1200

sel=(r0>=rmin)*(r0<=rmax)
popt, pcov = curve_fit(lin, r0[sel], w0[sel], sigma=ew0[sel])
w0small=w0-lin(r0, *popt)



# Write all results to disk
ascii.write([["fhg", "efhg", "Nh", "Ng", "M_GMC", "L_HII"], np.array([fhg, efhg, len(x2[sel2]), len(x1[sel1]), np.median(mgmc[sel1]), np.median(lha[sel2])])], './output/'+galaxy+'_fhg.txt', format='no_header', overwrite=True)
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
ax.axhline(y=0, linestyle=':', color='grey')
ax.axvline(x=rmax, linestyle=':', color='red')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'$\omega(r)$', fontsize=20)
ax.set_title(galaxy+" ; fhg="+"{:.2f}".format(fhg)+" ("+"{:.2f}".format(efhg)+")", fontsize=30)
ax.tick_params(labelsize=20)
#ax.set_yscale('log')
plt.savefig('./plots/'+galaxy+'_corr.png')
#plt.show()

fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(r0, w0, ew0, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0, w0, 'o', color='black', alpha=1.0)
#ax.set_xlim(0, 250)
ax.set_ylim(1e-3, 1e1)
ax.axhline(y=0, linestyle=':', color='grey')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'log($\omega(r)$)', fontsize=20)
ax.set_title(galaxy+" ; fhg="+"{:.2f}".format(fhg)+" ("+"{:.2f}".format(efhg)+")", fontsize=30)
ax.tick_params(labelsize=20)
ax.set_yscale('log')
plt.savefig('./plots/'+galaxy+'_corr_log.png')
#plt.show()

#Plot large scale fit and small scale cross-correlation function

fig, ax = plt.subplots(figsize=(12, 8))
ax.errorbar(r0, w0, ew0, fmt="o", color='grey', capsize=5, alpha=0.5)
ax.plot(r0, w0, 'o', color='black', alpha=1.0)
ax.plot(r0[sel], w0[sel], 'o', color='red', alpha=0.5)
ax.plot(r0[(r0<=rmax)], (lin(r0, *popt))[(r0<=rmax)], color='red')
ax.set_xlim(0, rmax)
ax.axhline(y=0, linestyle=':', color='grey')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'$\omega(r)$', fontsize=20)
ax.set_title(galaxy+" ; fhg="+"{:.2f}".format(fhg)+" ("+"{:.2f}".format(efhg)+")", fontsize=30)
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
ax.axhline(y=0, linestyle=':', color='grey')
plt.xlabel('r [pc]', fontsize=20)
plt.ylabel(r'$\omega(r)-\omega_{lin}$', fontsize=20)
plt.title(galaxy+" ; fhg="+"{:.2f}".format(fhg)+" ("+"{:.2f}".format(efhg)+")", fontsize=30)
ax.tick_params(labelsize=20)
#ax.set_yscale('log')
#ax.set_xscale('log')
plt.savefig('./plots/'+galaxy+'_corr_small.png')
#plt.show()




# ## Plot GMCs with HII regions and Random Catalog


fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(36, 12))

coflux=comap[0].data
#ax0.imshow(coflux, cmap='Blues', norm=LogNorm(), extent=[np.min(xmask), np.max(xmask), np.min(ymask), np.max(ymask)], alpha=0.5)
#ax0.imshow(hamap2, cmap='Reds', norm=LogNorm(), extent=[np.min(xmask), np.max(xmask), np.min(ymask), np.max(ymask)], alpha=0.5)

#ax0.contourf(xmask, ymask, coflux, cmap='Blues', norm=LogNorm(), alpha=0.3)
#ax0.contourf(xmask, ymask, hamap2, cmap='Reds', norm=LogNorm(), alpha=0.3)
ax0.contourf(xmask, ymask, np.log10(coflux), cmap='Blues', alpha=0.3, levels=50)
ax0.contourf(xmask, ymask, np.log10(hamap2), cmap='Reds', alpha=0.3, levels=50)
ax0.contour(xmask, ymask, comask2, colors='blue', alpha=0.5)
ax0.contour(xmask, ymask, hamask, colors='red', alpha=0.5)
ax0.contour(xmask, ymask, mask, colors='green', alpha=0.5)
#ax0.set_xlim(np.min(xmask[mask]), np.max(xmask[mask]))
#ax0.set_ylim(np.min(ymask[mask]), np.max(ymask[mask]))
ax0.set_aspect('equal')
ax0.set_xlabel('X [pc]', fontsize=20)
ax0.set_ylabel('Y [pc]', fontsize=20)
ax0.tick_params(labelsize=20)
custom_lines = [Line2D([0], [0], color='blue', lw=4), Line2D([0], [0], color='red', lw=4), Line2D([0], [0], color='green', lw=4)]
ax0.legend(custom_lines, ['ALMA', 'MUSE', "Mask"], fontsize=20)

ax1.plot(x1[sel1], y1[sel1], 'ob', alpha=0.3, label='GMCs')
ax1.plot(x2[sel2], y2[sel2], 'or', alpha=0.3, label='HII Regions')
ax1.plot(x1, y1, '.b', alpha=0.3)
ax1.plot(x2, y2, '.r', alpha=0.3)
ax1.set_xlim(np.min(xmask[mask]), np.max(xmask[mask]))
ax1.set_ylim(np.min(ymask[mask]), np.max(ymask[mask]))
ax1.contour(xmask, ymask, mask)
ax1.set_aspect('equal')
ax1.set_xlabel('X [pc]', fontsize=20)
ax1.set_ylabel('Y [pc]', fontsize=20)
ax1.set_title(galaxy, fontsize=30)
ax1.tick_params(labelsize=20)
ax1.legend(fontsize=20)

ax2.plot(x1[sel1], y1[sel1], 'ob', alpha=0.3, label='GMCs')
ax2.plot(x1, y1, '.b', alpha=0.3)
ax2.plot(xr, yr, '.g', alpha=0.1, label='Random Points')
ax2.set_xlim(np.min(xmask[mask]), np.max(xmask[mask]))
ax2.set_ylim(np.min(ymask[mask]), np.max(ymask[mask]))
ax2.contour(xmask, ymask, mask)
ax2.set_aspect('equal')
ax2.set_xlabel('X [pc]', fontsize=20)
ax2.set_ylabel('Y [pc]', fontsize=20)
ax2.tick_params(labelsize=20)
ax2.legend(fontsize=20)
plt.savefig('./plots/'+galaxy+'_xy.png')


#fig, ax = plt.subplots(figsize=(12, 12))
#ax.plot(x1[sel1], y1[sel1], 'ob', alpha=0.3, label='GMCs')
#ax.plot(x2[sel2], y2[sel2], 'or', alpha=0.3, label='HII Regions')
#ax.plot(x1, y1, '.b', alpha=0.3)
#ax.plot(x2, y2, '.r', alpha=0.3)
##ax.plot(xr, yr, '.g', alpha=0.1)
#ax.set_xlim(-1000, 1000)
#ax.set_ylim(-1000, 1000)
#ax.contour(xmask, ymask, mask)
#plt.xlabel('X [pc]', fontsize=20)
#plt.ylabel('Y [pc]', fontsize=20)
#plt.title(galaxy, fontsize=30)
#ax.tick_params(labelsize=20)
#ax.legend(fontsize=20, loc=1)
#plt.savefig('./plots/'+galaxy+'_xy_zoom.png')
##plt.show()
#
#fig, ax = plt.subplots(figsize=(12, 12))
#ax.plot(x1[sel1], y1[sel1], 'ob', alpha=0.3, label='GMCs')
##ax.plot(x2[sel2], y2[sel2], 'or', alpha=0.3)
#ax.plot(x1, y1, '.b', alpha=0.3)
##ax.plot(x2, y2, '.r', alpha=0.3)
#ax.plot(xr, yr, '.g', alpha=0.1, label='Random Points')
#ax.set_xlim(-1000, 1000)
#ax.set_ylim(-1000, 1000)
#ax.contour(xmask, ymask, mask)
#plt.xlabel('X [pc]', fontsize=20)
#plt.ylabel('Y [pc]', fontsize=20)
#plt.title(galaxy, fontsize=30)
#ax.tick_params(labelsize=20)
#ax.legend(fontsize=20)
#plt.savefig('./plots/'+galaxy+'_xy_rand_zoom.png')
##plt.show()





