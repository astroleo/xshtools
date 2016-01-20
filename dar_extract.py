from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import os
import glob
from astropy.modeling import models,fitting
from atmo_disp import *


def get_file(dpid):
	subdir=dpid.replace(":","_")+"_tpl"
	f=glob.glob(os.getenv("XDIRRED")+"/"+subdir+"/*_IFU_MERGE3D_DATA_OBJ_*.fits")
	if len(f) != 1:
		raise ValueError("There is no (or more than one) merged file in sub-directory " + subdir)
	
	return(f[0])

def fit_traces(flux_cube,wave):
	nslices=3

	##
	## number of wavelength pixels to sum for trace determination (has to be an even number)
#	nsum=100
	nsum=1000
	nwave=flux_cube.shape[0]
	nfit=int(np.floor(nwave/nsum))

	traces=np.zeros((nslices,nfit))

	for s in np.arange(nslices):
		sd=flux_cube[:,:,s]
		ny=sd.shape[1]

		for t in np.arange(nfit):
			g_init=models.Gaussian1D(amplitude=10000,mean=10,stddev=5)
			fit_g=fitting.LevMarLSQFitter()
			g=fit_g(g_init,np.arange(ny),np.sum(sd[t*nsum:(t+1)*nsum,:],axis=0))
			traces[s,t]=g.mean.value

	waves=wave[int(nsum/2)+nsum*np.arange(nfit)]
	return(waves,traces)

##
## purpose: give the DPID of a NIR observation, get maximum position, compute
##    pixel position in UVB and VIS arm taking into account differential
##    atmospheric refraction
##
##    if also DPID_VIS and DPID_UVB are given, compute traces for all three arms
##       and overplot computed position over actual position, saves plot to file

def dar_position(DPID_NIR,DPID_VIS=False,DPID_UVB=False):
	
	if (DPID_VIS != False) and (DPID_UVB == False):
		raise ValueError("If DPID_VIS is given, you must also give DPID_UVB")
	if (DPID_UVB != False) and (DPID_VIS == False):
		raise ValueError("If DPID_UVB is given, you must also give DPID_VIS")
	
	if DPID_VIS != False:
		QCplot=True
	else:
		QCplot=False

	fNIR=get_file(DPID_NIR)
	if QCplot==True:
		fUVB=get_file(DPID_UVB)
		fVIS=get_file(DPID_VIS)

	##
	## get relevant data from files
	hdu_NIR=fits.open(fNIR)
	flux_NIR=hdu_NIR[0].data
	hdr_NIR=hdu_NIR[0].header
	wave_NIR=(np.arange(hdr_NIR['NAXIS3'])+1-hdr_NIR['CRPIX3'])*hdr_NIR['CDELT3']+hdr_NIR['CRVAL3']
	CRPIX2_NIR=hdr_NIR['CRPIX2']
	CDELT2_NIR=hdr_NIR['CDELT2']

	if QCplot==True:
		hdu_UVB=fits.open(fUVB)
		flux_UVB=hdu_UVB[0].data
		hdr_UVB=hdu_UVB[0].header

		hdu_VIS=fits.open(fVIS)
		flux_VIS=hdu_VIS[0].data
		hdr_VIS=hdu_VIS[0].header

	##
	## hard-coded for the moment, so that I don't have to find and read UVB and VIS files
	## to be replace by private class variables when all of this becomes a class some time...
	naxis3_UVB=14554
	crpix3_UVB=0
	cdelt3_UVB=0.0199999995529652
	crval3_UVB=298.903594970703
	wave_UVB=(np.arange(naxis3_UVB)+1-crpix3_UVB)*cdelt3_UVB+crval3_UVB

	naxis3_VIS=24317
	crpix3_VIS=0
	cdelt3_VIS=0.0199999995529652
	crval3_VIS=533.657836914062
	wave_VIS=(np.arange(naxis3_VIS)+1-crpix3_VIS)*cdelt3_VIS+crval3_VIS

	CRPIX2_VIS=13
	CDELT2_VIS=4.44444434510337e-05

	##
	## determine center position in NIR arm
	## use median of J and H slices (since K band has weird bump problem)
	J=(wave_NIR>1150) & (wave_NIR<1350)
	H=(wave_NIR>1509) & (wave_NIR<1799)
	profile=np.sum(flux_NIR[J,:,:],axis=(0,2)) + np.sum(flux_NIR[H,:,:],axis=(0,2))
	xx=np.arange(len(profile))

	g_init=models.Gaussian1D(amplitude=100000,mean=10,stddev=5)
	fit_g=fitting.LevMarLSQFitter()
	g=fit_g(g_init,xx,profile)
	y_NIR=g.mean.value

	##
	## convert pixel position from NIR to red end of VIS (where differential
	##    refraction essentially does not play a role), CRPIX2 and CDELT2 are
	##    different between UVB/VIS and NIR arms.
	y_VIS = (y_NIR-CRPIX2_NIR) + CRPIX2_VIS
	px_scale_arcsec_UVBVIS = 3600 * CDELT2_VIS
	##
	## determine differential atmospheric refraction (DAR) correction from NIR to UVB
	wave_ref=1000
	waves_all = 300 + np.arange(2000)
	waves_all_mu = waves_all/1000
	zeta_deg = 90-hdr_NIR['HIERARCH ESO TEL ALT']
	##
	## dd_pixel input is in micron
	dd_pixel_VIS = y_VIS + dDAR(zeta_deg, wave_VIS/1000, wave_ref/1000)/px_scale_arcsec_UVBVIS
	##
	## use y_VIS offset because reference position is same as for dd_pixel_VIS above
	dd_pixel_UVB = y_VIS + dDAR(zeta_deg, wave_UVB/1000, wave_ref/1000)/px_scale_arcsec_UVBVIS

	if QCplot==True:
		twaves_NIR,traces_NIR=fit_traces(flux_NIR,wave_NIR)
		twaves_VIS,traces_VIS=fit_traces(flux_VIS,wave_VIS)
		twaves_UVB,traces_UVB=fit_traces(flux_UVB,wave_UVB)

		plt.plot(twaves_NIR,traces_NIR[0,:],'b.')
		plt.plot(twaves_NIR,traces_NIR[1,:],'g.')
		plt.plot(twaves_NIR,traces_NIR[2,:],'r.')

		plt.plot(twaves_VIS,traces_VIS[0,:],'b.')
		plt.plot(twaves_VIS,traces_VIS[1,:],'g.')
		plt.plot(twaves_VIS,traces_VIS[2,:],'r.')

		plt.plot(twaves_UVB,traces_UVB[0,:],'b.')
		plt.plot(twaves_UVB,traces_UVB[1,:],'g.')
		plt.plot(twaves_UVB,traces_UVB[2,:],'r.')

		plt.plot(wave_VIS,dd_pixel_VIS)
		plt.plot(wave_UVB,dd_pixel_UVB)
		plt.plot([1000,2500],[y_NIR,y_NIR])
		plt.xlim([200,2500])
		plt.ylim([0,30])
		fplot=os.getenv("XSHOOANALYSIS")+"/QC/plot_DAR/"+DPID_NIR.replace(":","_")+".png"
		plt.savefig(fplot)
		plt.clf()

	return(dd_pixel_UVB, dd_pixel_VIS, y_NIR)


#UVB="XSHOO.2014-02-25T08:57:49.298"
#VIS="XSHOO.2014-02-25T08:57:54.450"
#NIR="XSHOO.2014-02-25T08:57:57.817"

## altitude: 34deg
#NIR="XSHOO.2015-11-19T04:05:32.143"
#UVB="XSHOO.2015-11-19T04:05:23.614"
#VIS="XSHOO.2015-11-19T04:05:28.804"

#wave_UVB,wave_VIS,dd_pixel_UVB,dd_pixel_VIS = dar_position(NIR,DPID_VIS=VIS,DPID_UVB=UVB)
#wave_UVB,wave_VIS,dd_pixel_UVB,dd_pixel_VIS = dar_position(NIR)
