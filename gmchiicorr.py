import time
import numpy as np
from scipy.spatial.distance import cdist
from fast_histogram import histogram1d
import astropy.io.fits as fits
from scipy.interpolate import NearestNDInterpolator
import random



# ## Routine to get basic positional parameters of a galaxy

def getparams(name):
    sample=fits.open('./catalogs/phangs_sample_table_v1p4.fits')
    names=sample[1].data['NAME']
    sel=sample[1].data['NAME']==name
    data=sample[1].data[sel]
    
    ra0=data['ORIENT_RA'][0]
    dec0=data['ORIENT_DEC'][0]
    D=data['DIST'][0]
    PA=data['ORIENT_POSANG'][0]
    inc=data['ORIENT_INCL'][0]
    
    return (ra0, dec0, D, PA, inc)


# Fucntion that converts RA,DEC into x,y in pc
# ra, dec is input
# ra0, dec0 is center of the galaxy
# D is distance in Mpc
# PA is position angle
# i is inclination
# return x, y, r, theta (i.e. positions on disk in cartesian and polar coordinates)
def radec2xy(ra, dec, ra0, dec0, D, PA, inc):
    dra=(ra0-ra)*3600*np.cos(np.pi/180*dec) # relative ra from center [arcsec], *-1 so E-->right
    ddec=(dec-dec0)*3600.                  # relative dec from center [arcsec]
    dx1=dra*np.sin(np.pi/180*PA)+ddec*np.cos(np.pi/180*PA)
    dy1=-1.*dra*np.cos(np.pi/180*PA)+ddec*np.sin(np.pi/180*PA)
    rsec=np.sqrt((dy1/np.cos(np.pi/180*inc))**2+dx1**2)
    theta=np.arctan(dy1/dx1)*180/np.pi
    theta[dx1<=0]=theta[dx1<=0]+180.
    theta[(dx1>=0)*(dy1<0)]=theta[(dx1>=0)*(dy1<0)]+360.
    pc2sec=D*1e6*np.pi/180/3600.  # pc/arcsec at distance D

    return (dx1*pc2sec,dy1/np.cos(np.pi/180*inc)*pc2sec,rsec*pc2sec,theta, dra, ddec)

# ## Function that makes random catalog of N objects inside mask area

def rancat(xmask, ymask, mask, Nr=10000, f=None):
    xr=np.random.uniform(np.min(xmask), np.max(xmask), 1000)
    yr=np.random.uniform(np.min(ymask), np.max(ymask), 1000)
    if f == None:
        aux=np.stack((xmask.flatten(), ymask.flatten()), axis=-1)
        f=NearestNDInterpolator(aux, mask.flatten())
    selr=f(xr, yr)
    frac=len(selr[selr])/len(selr)
    xr=np.random.uniform(np.min(xmask), np.max(xmask), int(Nr/frac))
    yr=np.random.uniform(np.min(ymask), np.max(ymask), int(Nr/frac))
    selr=f(xr, yr)
    return(xr[selr], yr[selr])


# Function to calculate binomial errors for ratios (needs ratio<1)
import scipy.stats.distributions as dist
# k = success count
# n = sample size
# c = confidence level (c=0.683 == 1 sigma for Normal distribution)
def ebratio(k, n, c=0.683):
    p_lower=dist.beta.ppf((1-c)/2.,k+1,n-k+1)
    p_upper=dist.beta.ppf(1-(1-c)/2.,k+1,n-k+1)
    return(k/n, k/n-p_lower, p_upper-k/n)

# Function to calculate 1-sigma gaussian errror for ratios (needs k and n >> 1) # PREFERRED
# using "Uncorrelated noncentral normal ratio" formula
def egratio(k, n):
    sigk=np.sqrt(k)
    sign=np.sqrt(n)
    ratio=k/n
    sigratio=(k/n)*np.sqrt((sigk/k)**2+(sign/n)**2)
    return(ratio, sigratio)

# Function to calculate 1-sigma gaussian errror for ratios (needs k and n >> 1)
# using "Cameron 2011 from Elmegreen 1990" formula
def enratio(k, n):
    p=k/n
    q=1-p
    sigp=p*np.sqrt(p*q/n)
    return(p, sigp)

# ## Function that calculates the 2-point cross-correlation function
# "bins" must be centers of bins and regularly spaced if provided by user

def w(x1, y1, x2, y2, xr, yr, rmin=0, rmax=5000, dr=25, bins=None):

    c1=np.transpose(np.array([x1,y1]))
    c2=np.transpose(np.array([x2,y2]))
    cr=np.transpose(np.array([xr,yr]))

    ddarr=cdist(c1, c2).ravel()  
    drarr=cdist(c1, cr).ravel()

#    ddarr=fastdist.matrix_to_matrix_distance(c1, c2, fastdist.euclidean, "euclidean")
#    drarr=fastdist.matrix_to_matrix_distance(c1, cr, fastdist.euclidean, "euclidean")

    # Count pairs in distance bins       
    N2=len(x2)
    N1=len(x1)
    Nr=len(xr)

    # if no "bins" array is provided use rmin, rmax, and dr to compute it, if "bins" is provided compute dr
    if bins is not None:
        dr=bins[1]-bins[0]
    else:
        bins=np.arange(rmin+dr/2, rmax, dr) # centers of bins for output

        
#    dd0, dd0bins = np.histogram(ddarr, bins=np.arange(rmin, rmax+dr, dr)) #here bins are bin edges
#    dr0, dr0bins = np.histogram(drarr, bins=np.arange(rmin, rmax+dr, dr))

    dd0 = histogram1d(ddarr, bins=len(bins), range=(bins[0]-dr/2,bins[-1]+dr/2))  #here range are bin edges
    dr0 = histogram1d(drarr, bins=len(bins), range=(bins[0]-dr/2,bins[-1]+dr/2))
    
    # Normalize pair counts and compute cross-correlation function
    if (N1!=0)*(N2!=0)*(Nr!=0):
        dd=dd0/N1/N2
        dr=dr0/N1/Nr
        omega=dd/dr-1
    
        edd=np.sqrt(dd0)/N1/N2
        edr=np.sqrt(dr0)/N1/Nr

        eomega=np.sqrt((edd/dr)**2+(dd*edr/dr**2)**2)
    else:
        omega=np.repeat(np.nan, len(bins))
        eomega=np.repeat(np.nan, len(bins))
        
    return (bins, omega, eomega)


# ## Functions that draws N model GMCs and returns properties

# This version just uses the observed GMC coordinates and assumes constant parameters for all clouds

def drawgmc_xy(x,y,rc=25,tc=30,ts=10,tfb=5,Ng=1,voff=10):

    ngmc=len(x) # number of GMCs
    
    xgmc=x # gmc x coordinates [pc]
    ygmc=y # gmc y coordinates [pc]
    
    rc=np.repeat(rc,ngmc) # gmc radius [pc]
    tc=np.repeat(tc,ngmc) # gmc lifetime [Myr]
    ts=np.repeat(ts,ngmc) # stellar SF tracer lifetime [Myr]
    tfb=np.repeat(tfb,ngmc) # stellar SF tracer emergence time [Myr]
    Ng=np.repeat(Ng,ngmc).astype(int) # number of massive SF episodes during gmc lifetime
    voff=np.repeat(voff,ngmc) # stellar SF tracer centorid to GMC centroid offset speed [km/s]

    tobs=np.random.uniform(0, tc-tfb+ts, ngmc) # random observation time for this cloud (between zero and tc+ts)
    fgmc=np.repeat(True, ngmc) # GMC visibility flag
    fgmc[tobs > tc]=False # visibility flag is False if GMC has faded

    return (xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc)
    

# This version samples GMC coordinates with a typical separation scale "l" within a dbox**2 kpc**2 box, and assumes constant parameters for all clouds
# Also returns random catalog coordinates over the same area and Frand*Ngmc points
def drawgmc_l(dbox=2000,l=200,rc=25,tc=30,ts=10,tfb=5,Ng=1,voff=10, frand=10):

    area=dbox**2
    ngmc=int(area/l**2) # number of GMCs
    
    xgmc=np.random.uniform(-0.5*dbox, 0.5*dbox, ngmc) # gmc x coordinates [pc]
    ygmc=np.random.uniform(-0.5*dbox, 0.5*dbox, ngmc) # gmc y coordinates [pc]

    xr=np.random.uniform(-0.5*dbox, 0.5*dbox, frand*ngmc) # gmc x coordinates [pc]
    yr=np.random.uniform(-0.5*dbox, 0.5*dbox, frand*ngmc) # gmc y coordinates [pc]

    rc=np.repeat(rc,ngmc) # gmc radius [pc]
    tc=np.repeat(tc,ngmc) # gmc lifetime [Myr]
    ts=np.repeat(ts,ngmc) # stellar SF tracer lifetime [Myr]
    tfb=np.repeat(tfb,ngmc) # stellar SF tracer emergence time [Myr]
    Ng=np.repeat(Ng,ngmc).astype(int) # number of massive SF episodes during gmc lifetime
    voff=np.repeat(voff,ngmc) # stellar SF tracer centorid to GMC centroid offset speed [km/s]
    tobs=np.random.uniform(0, tc-tfb+ts, ngmc) # random observation time for this cloud (between zero and tc+ts)
    fgmc=np.repeat(True, ngmc) # GMC visibility flag
    fgmc[tobs > tc]=False # visibility flag is False if GMC has faded

    #print("Generating GMC and random coordinates:", ngmc, frand*ngmc)
    return (xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc, xr, yr)
    


# ## Function that draws N*Ng HII regions given an ensemble of GMCs
#
#def drawhii(xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc):
#   
#   ngmc=len(xgmc) # number of GMCs
#   nhii=np.sum(Ng.astype(int)) # number of HII regions
#   voffaux=voff*3.2e-14*(1e6*365*24*3600) # velocity in pc/Myr
#
#   xhii=np.zeros(nhii) # HII region initial coordinates
#   yhii=np.zeros(nhii)
#   t0=np.zeros(nhii) # HII region formation time
#   fhii=np.repeat(False,nhii) # HII region visibility flag
#   indgmc=np.zeros(nhii) # gmc list index
#
#   k=0 # counter
#   for i in range(ngmc):  # drawing HII regions for each GMC independently
#       for j in range(Ng[i]):
#
#           # initial position
#           rad0=rc[i]*np.sqrt(random.uniform(0,1)) # uniform across circular area
#           #rad0=rc[i]*np.random.uniform(0,1) # uniform in radius (i.e. as r**-2)
#           theta0=random.uniform(0, 2*np.pi)
#           x0=xgmc[i]+rad0*np.cos(theta0)
#           y0=ygmc[i]+rad0*np.sin(theta0)
#
#           #formation time (assuming stars form no later than tc-tfb)
#           t0[k]=random.uniform(0, tc[i]-tfb[i])    
#           
#           #offset direction and final position at tobs
#           phi=random.uniform(0, 2*np.pi) # random velocity angle on plane of the galaxy
#           xhii[k]=x0+voffaux[i]*(tobs[i]-t0[k])*np.cos(phi)
#           yhii[k]=y0+voffaux[i]*(tobs[i]-t0[k])*np.sin(phi)
#
#           #visibility flag is True if cloud has already formed and emerged, and has not yet faded
#           if (t0[k]+tfb[i]<tobs[i])*(tobs[i]<t0[k]+ts[i]):
#               fhii[k]=True
#
#               #print(t0[k], tfb[i], ts[i], tobs[i], t0[k]+tfb[i], t0[k]+ts[i], fhii[k])
#           
#           #GMC list index
#           indgmc[i]=i
#           
#           k=k+1
#       
#   
#   return (xhii,yhii,fhii)
    
# # ## Function that draws N*Ng HII regions given an ensemble of GMCs # Faster version!
#
#def drawhii(xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc):
#    
#    ngmc=len(xgmc) # number of GMCs
#    nhii=np.sum(Ng.astype(int)) # number of HII regions
#    voffaux=voff*3.2e-14*(1e6*365*24*3600) # velocity in pc/Myr
#
#    rc0=np.zeros(nhii) # 
#    xgmc0=np.zeros(nhii)
#    ygmc0=np.zeros(nhii)
#    voffaux0=np.zeros(nhii) 
#    tobs0=np.zeros(nhii) 
#    ts0=np.zeros(nhii) 
#    tfb0=np.zeros(nhii) 
#
#    xhii=np.zeros(nhii) # HII region initial coordinates
#    yhii=np.zeros(nhii)
#    t0=np.zeros(nhii) # HII region formation time
#    fhii=np.repeat(False,nhii) # HII region visibility flag
#    indgmc=np.zeros(nhii) # gmc list index
#
#    k=0 # counter
#    for i in range(ngmc):  # drawing HII regions for each GMC independently
#        for j in range(Ng[i]):
#            rc0[k]=rc[i]        # GMC radius for this HII region
#            xgmc0[k]=xgmc[i]    # GMC x for this HII region
#            ygmc0[k]=ygmc[i]     # GMC y for this HII region
#            voffaux0[k]=voffaux[i]
#            tobs0[k]=tobs[i]
#            ts0[k]=ts[i]
#            tfb0[k]=tfb[i]
#            #formation time (assuming stars form no later than tc-tfb)
#            t0[k]=random.uniform(0, tc[i]-tfb[i])
#            indgmc[i]=i         #GMC list index
#            k=k+1
#
#    # initial position
#    rad0=rc0*np.sqrt(np.random.uniform(0,1,nhii)) # uniform across circular area
##   rad0=rc0*np.random.uniform(0,1, nhii) # uniform in radius (i.e. as r**-2)
#    theta0=np.random.uniform(0, 2*np.pi,nhii)
#    x0=xgmc0+rad0*np.cos(theta0)
#    y0=ygmc0+rad0*np.sin(theta0)
#     
#    #offset direction and final position at tobs
#    phi=np.random.uniform(0, 2*np.pi,nhii) # random velocity angle on plane of the galaxy
#    xhii=x0+voffaux0*(tobs0-t0)*np.cos(phi)
#    yhii=y0+voffaux0*(tobs0-t0)*np.sin(phi)
#
#     #visibility flag is True if cloud has already formed and emerged, and has not yet faded
#    fhii[(t0+tfb0<tobs0)*(tobs0<t0+ts0)]=True
#
#
#    
#    return (xhii,yhii,fhii)
    
# ## Function that draws N*Ng HII regions given an ensemble of GMCs # FAstest!! but not working well

def drawhii(xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc):
    
    ngmc=len(xgmc) # number of GMCs
    nhii=np.sum(Ng.astype(int)) # number of HII regions
    voffaux=voff*3.2e-14*(1e6*365*24*3600) # velocity in pc/Myr

    xgmc0=np.zeros(nhii)
    ygmc0=np.zeros(nhii)
    tobs0=np.zeros(nhii) 

    k=0 # counter
    for i in range(ngmc):  # drawing HII regions for each GMC independently
        xgmc0[k:k+Ng[i]-1]=xgmc[i]    # GMC x for this HII region
        ygmc0[k:k+Ng[i]-1]=ygmc[i]     # GMC y for this HII region
        tobs0[k:k+Ng[i]-1]=tobs[i]     # GMC y for this HII region
        k=k+Ng[i]

    
    #formation time (assuming stars form no later than tc-tfb)
    t0=np.random.uniform(0, tc[0]-tfb[0], nhii)

    # initial position
    rad0=rc[0]*np.sqrt(np.random.uniform(0,1,nhii)) # uniform across circular area
#   rad0=rc0*np.random.uniform(0,1, nhii) # uniform in radius (i.e. as r**-2)
    theta0=np.random.uniform(0, 2*np.pi,nhii)
    x0=xgmc0+rad0*np.cos(theta0)
    y0=ygmc0+rad0*np.sin(theta0)
     
    #offset direction and final position at tobs
    phi=np.random.uniform(0, 2*np.pi,nhii) # random velocity angle on plane of the galaxy
    xhii=x0+voffaux[0]*(tobs0-t0)*np.cos(phi)
    yhii=y0+voffaux[0]*(tobs0-t0)*np.sin(phi)

    
     #visibility flag is True if cloud has already formed and emerged, and has not yet faded
    fhii=np.repeat(False,nhii) # HII region visibility flag
    fhii[(t0+tfb[0]<tobs0)*(tobs0<t0+ts[0])]=True


    
    return (xhii,yhii,fhii)
         

# ## Linear Model for fitting large scale correlation function

def lin(x, a, b):
    return a+b*x


# ## Function that Evaluates Cross Correlation Function for Model Parameters

def eval_w(l0, rc0, tc0, ts0, tfb0, Ng0, voff0, rmin=0, rmax=500, dr=25, bins=False, Nsamples=500):
    
    t0=time.time()
    
    print("Evaluating Model:")
    print("l=", l0)
    print("rc=", rc0)
    print("tc=", tc0)
    print("ts=", ts0)
    print("tfb=", tfb0)
    print("Ng=", Ng0)
    print("voff=", voff0)

    
    # Run w() one time to get bins
    xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc, xr, yr = drawgmc_l(dbox=2000, l=l0, rc=rc0, tc=tc0, ts=ts0, tfb=tfb0, Ng=Ng0, voff=voff0, frand=10)
    xhii, yhii, fhii = drawhii(xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc)
    print("Number of GMCs and HII Regions", len(xgmc[fgmc]), len(xhii[fhii]))
    print("fhg=", len(xhii[fhii])/len(xgmc[fgmc]))
    r0, w0, ew0 = w(xgmc[fgmc], ygmc[fgmc], xhii[fhii], yhii[fhii], xr, yr, rmin=rmin, rmax=rmax, dr=dr, bins=bins)   
    w0arr=np.zeros((len(w0),Nsamples))
    ew0arr=np.zeros((len(ew0),Nsamples))
    fhgarr=np.zeros(Nsamples)  # ratio of number of observed hii regions to gmcs
    
   
    # Run w() Nsample times and average
    for i in range(Nsamples):
    
        xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc, xr, yr = drawgmc_l(dbox=2000, l=l0, rc=rc0, tc=tc0, ts=ts0, tfb=tfb0, Ng=Ng0, voff=voff0, frand=10)
        xhii, yhii, fhii = drawhii(xgmc, ygmc, rc, tc, ts, tfb, Ng, voff, tobs, fgmc)
        if (len(xgmc[fgmc])!=0)*(len(xhii[fhii])!=0):
            r0, w0arr[:,i], ew0arr[:,i] = w(xgmc[fgmc], ygmc[fgmc], xhii[fhii], yhii[fhii], xr, yr, rmin=rmin, rmax=rmax, dr=dr, bins=bins)
            fhgarr[i]=len(xhii[fhii])/len(xgmc[fgmc])
        else:
            w0arr[:,i]=np.repeat(np.nan,len(r0))
            ew0arr[:,i]=np.repeat(np.nan,len(r0))
            fhgarr[i]=np.nan
            
    w0arr[~np.isfinite(w0arr)]=np.nan
    ew0arr[~np.isfinite(ew0arr)]=np.nan    
    w0=np.nanmean(w0arr, axis=1)
    ew0=np.nanmean(ew0arr, axis=1)/np.sqrt(Nsamples)
    fhg0=np.nanmean(fhgarr)

    print(w0)  
    print("Model Evaluation Run Time [s] =", time.time()-t0)
  
    return (r0, w0, ew0, fhg0)
