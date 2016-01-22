##
## X-SHOOTER data post processing
##
## (1) D.A.R. correction and 1D spectral extraction of spec + err

## (2?) telluric and flux calibration
## (3?) stitch spectrum together / find fudge factors between arms
##

#
## required environment variables
##
## XDIR -- root directory for QC plots and output (calibrated) files
## XDIRRED -- root directory where the pipeline-reduced data are

import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
from matplotlib import pyplot as plt
from dar_extract import dar_position
import sqlite3

import os
import pdb

	
##
## FUNCTION qual_interpret_interpolate
##
## PURPOSE
##    - interpret QUAL extension in X-SHOOTER reduced data
##
## INPUT
##     flux_cube, noise_cube
##     qual_cube: cube with bad pixel codes
##     ix: aperture mask
##
def qual_interpret_interpolate(flux_cube, noise_cube, qual_cube, ix):
	## code_bad: code above which a pixel is interpreted as bad
	code_bad = 1

	##
	## bad flags within aperture (bad_aper) and within aperture + 1 pixel (bad_aper2) -- for interpolation purposes
	bad_aper = ix*qual_cube >= code_bad
	
	qual=np.zeros(flux_cube.shape[0])
	corr=np.zeros(flux_cube.shape[0])

	ncorr=0
	nbad=0
	ngood=0
	
	for s in np.arange(flux_cube.shape[0]):
		nbad_in_aper = np.sum(bad_aper[s,:,:])
		
		if nbad_in_aper == 0:
			qual[s] = 1
			corr[s] = 0
			ngood+=1
		else:
			qual[s] = 0
			corr[s] = 0
			nbad+=1


	print("bad: {0}, good: {1} pixels.".format(nbad,ngood))
	flux = np.sum(flux_cube*ix, axis=(1,2))
	##
	## "the variance of a sum of uncorrelated events is the sum of the variances"
	noise = np.sqrt(np.sum((noise_cube*ix)**2, axis=(1,2)))
		
	return(flux,noise,qual,corr)


##
## FUNCTION flatten_spectrum
## 
## PURPOSE
##    take data cube, extract 1D spectrum from it
##
def flatten_spectrum(f,arm,dd_pixel,outfile=None):
	hdu=fits.open(f)
	##
	## wavelength vector is the same for all three extensions, 
	## i.e. we can get it from main header
	##
	hdr=hdu[0].header
	wave = (hdr['CRVAL3'] - hdr['CRPIX3'] * hdr['CDELT3'] + 
		hdr['CDELT3'] * (1+np.arange(hdr['NAXIS3'])))
	##
	flux_cube=hdu[0].data
	noise_cube=hdu[1].data
	qual_cube=hdu[2].data

	scale_y_arcsec=3600*float(hdu[0].header['CDELT2'])
	aperture_arcsec=1.0

	r_px_y = np.ceil(aperture_arcsec/(2*scale_y_arcsec))

	##
	## dd_pixel is either an array of same length as wave (check) or a single value (y_NIR from dar_position)
	if dd_pixel.size > 1:
		z,y,x=np.mgrid[:flux_cube.shape[0],:flux_cube.shape[1],:flux_cube.shape[2]]
		ix  = np.abs(y-dd_pixel[:,None,None]) < r_px_y
	else:
		z,y,x=np.mgrid[:flux_cube.shape[0],:flux_cube.shape[1],:flux_cube.shape[2]]
		ix  = np.abs(y-dd_pixel) < r_px_y

	##
	##
	flux,noise,qual,corr = qual_interpret_interpolate(flux_cube, noise_cube, qual_cube, ix)
	
	## plot (flux vs. wavelength + quality)
	norm=np.median(flux)
	plt.plot(wave,flux/norm)
	plt.plot(wave,noise/norm,'g-')
	plt.plot(wave[qual==0],flux[qual==0]/norm,'rx')
	plt.plot(wave[corr==1],flux[corr==1]/norm,'g.',ms=10)
	plt.xlabel("Wavelength")
	plt.ylabel("Flux (normalized)")
	if arm == "UVB":
		plt.ylim([0,3])
	else:
		plt.ylim([0,2])

	if outfile:		
		plt.savefig(outfile+".pdf")
		plt.close()

		prihdu=fits.PrimaryHDU(header=hdr)
		col1=fits.Column(name="WAVE", format='1E', array=wave)
		col2=fits.Column(name="FLUX", format='1E', array=flux)
		col3=fits.Column(name="NOISE", format='1E', array=noise)
		col4=fits.Column(name="QUAL", format='1J', array=qual)
		cols=fits.ColDefs([col1,col2,col3,col4])
		tbhdu=fits.BinTableHDU.from_columns(cols)
		tbhdulist = fits.HDUList([prihdu, tbhdu])
		tbhdulist.writeto(outfile,clobber=True)
		
		if arm == "VIS":
			## plot closeup in CaT region
			wregion=[850,880]		
		elif arm == "UVB":
			## plot closeup of H-beta 486 region
			wregion=[470,500]
		elif arm == "NIR":
			## plot closeup of Pa-beta 1282 region
			wregion=[1250,1350]
		
#		pdb.set_trace()
		norm=np.median(flux[(wave > wregion[0]) & (wave < wregion[1])])

		plt.plot(wave,flux/norm)
		plt.plot(wave,noise/norm,'g-')
		plt.plot(wave[qual==0],flux[qual==0]/norm,'rx',ms=2)
		plt.plot(wave[corr==1],flux[corr==1]/norm,'g.',ms=5)
		plt.xlabel("Wavelength")
		plt.ylabel("Flux (normalized)")
		plt.ylim([0,1.5])
		plt.xlim(wregion)
		plt.savefig(outfile+"_"+arm+"_detail.pdf")
		plt.close()
	else:
		plt.show()
		
	return(wave,flux,noise,qual,corr)

##
## helper function to inspect flux and mask cube and spectrum around chosen slice
##
## ...to be turned into an interactive widget at some point
##
## keyword wave: find slice closest to this wavelength
##
def inspect_cube(f,s,w=False):
	plt.close()
	hdu=fits.open(f)
	##
	## wavelength vector is the same for all three extensions, 
	## i.e. we can get it from main header
	##
	hdr=hdu[0].header
	wave = (hdr['CRVAL3'] - hdr['CRPIX3'] * hdr['CDELT3'] + 
		hdr['CDELT3'] * (1+np.arange(hdr['NAXIS3'])))

	##
	## if keyword wave is given try to find slice closest to wavelength
	if w == True:
		dl=hdr['CDELT3']
		ix=np.where((wave > (s - 0.5*dl)) & (wave < (s + 0.5*dl)))
		ss=ix[0][0]
		print("Found wavelength {0} at slice {1} (exact wavelength: {2})".format(s,ss,wave[ss]))
		s=ss	

	##
	flux_cube=hdu[0].data
#	noise_cube=hdu[1].data
	qual_cube=hdu[2].data
	##
	##
	plt.subplot2grid((1,7), (0,0), rowspan=1, colspan=1)
	plt.imshow(flux_cube[s,:,:],origin="lower",interpolation="nearest")
#	plt.colorbar()

	plt.subplot2grid((1,7), (0,1), rowspan=1, colspan=1)
	plt.imshow(qual_cube[s,:,:],origin="lower",vmin=0,vmax=32,interpolation="nearest")
	plt.colorbar()

	plt.subplot2grid((1,7), (0,2), rowspan=1, colspan=5)
	s1=np.max([0,s-200])
	s2=np.min([hdr['NAXIS3'],s+200])
	spec=np.sum(flux_cube,axis=(1,2))

##
## plot by wavelength
#	plt.plot(wave,spec/np.median(spec))
#	plt.plot([wave[s],wave[s]], [0,1], 'k.-.')
#	plt.xlim([wave[s1],wave[s2]])

##
## plot by slice
	slices=np.arange(len(spec))
	plt.plot(slices,spec/np.median(spec[s1:s2]))
	plt.plot([s,s], [0,2], 'k.-.')
	plt.xlim([s1,s2])

	plt.ylim([0.7,1.3])
	plt.tight_layout()
	plt.show()

##
## FUNCTION flatten_ob
##
## PURPOSE
##    prepare all data needed to calibrate one science observation
##
## INPUT
##    ob_name
##
def flatten_ob(ob_name):
	arms=["NIR","VIS","UVB"]
	for arm in arms:
		dataset_definition = os.getenv("XDIR")+'/dataset_definition/'+ob_name+'_'+arm+'.txt'
		a=ascii.read(dataset_definition,data_start=0)
		
		##
		## UVB does not need telluric correction
		if arm=="UVB":
			dprlist=["SCI","FLUX"]
		else:
			dprlist=["SCI","TELL","FLUX"]
		
		for row in a:
			ob=row[0]
			dpid=row[1]
			dpr=dprlist[row.index]
			dir_cube=os.getenv("XDIRRED")+"/"+dpid.replace(":","_")+"_tpl/"
			f_cube=dir_cube+ob+"_"+dpr+"_IFU_MERGE3D_DATA_OBJ_"+arm+".fits"

			if not os.path.isfile(f_cube):
				print("Cube file does not exist at", f_cube)
				continue

			dir_out=os.getenv("XDIR")+"/spectra/"+ob_name+"/"
			if not os.path.isdir(dir_out):
				os.mkdir(dir_out)
		
			f_out=dir_out+dpr+"_"+arm+"_spec.fits"
			if os.path.isfile(f_out):
				print("Outfile (spectrum) exists:",f_out)
				continue

			##
			## do centroiding on near-IR arm and compute atmospheric dispersion correction from that position
			if arm=="NIR" and dpr=="SCI":
				dd_pixel_UVB_SCI,dd_pixel_VIS_SCI,y_NIR = dar_position(dpid)
				flatten_spectrum(f_cube,arm,y_NIR,outfile=f_out)
			if arm=="NIR" and dpr=="TELL":
				dd_pixel_UVB_TELL,dd_pixel_VIS_TELL,y_NIR = dar_position(dpid)
				flatten_spectrum(f_cube,arm,y_NIR,outfile=f_out)
			if arm=="NIR" and dpr=="FLUX":
				dd_pixel_UVB_FLUX,dd_pixel_VIS_FLUX,y_NIR_FLUX = dar_position(dpid)
				flatten_spectrum(f_cube,arm,y_NIR,outfile=f_out)
			
			if arm=="UVB" and dpr=="SCI":
				flatten_spectrum(f_cube,arm,dd_pixel_UVB_SCI,outfile=f_out)
			if arm=="UVB" and dpr=="TELL":
				flatten_spectrum(f_cube,arm,dd_pixel_UVB_TELL,outfile=f_out)
			if arm=="UVB" and dpr=="FLUX":
				flatten_spectrum(f_cube,arm,dd_pixel_UVB_FLUX,outfile=f_out)

			if arm=="VIS" and dpr=="SCI":
				flatten_spectrum(f_cube,arm,dd_pixel_VIS_SCI,outfile=f_out)
			if arm=="VIS" and dpr=="TELL":
				flatten_spectrum(f_cube,arm,dd_pixel_VIS_TELL,outfile=f_out)
			if arm=="VIS" and dpr=="FLUX":
				flatten_spectrum(f_cube,arm,dd_pixel_VIS_FLUX,outfile=f_out)




################################################################################
#################################### "MAIN" ####################################
################################################################################

#dir="/Volumes/astrodata2/XSHOODATA-Reduced/reflex_end_products/2015-11-12_pipeline_2-6-8/"
#files=ascii.read("/Users/leo/Desktop/files.txt")
#for f in files:
#	ff=dir+f[0]
#	print(ff)
#	w,f,n,q,c = flatten_spectrum(ff)
