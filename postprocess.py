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

	nbad=0
	ngood=0
	
	for s in np.arange(flux_cube.shape[0]):
		nbad_in_aper = np.sum(bad_aper[s,:,:])
		
		if nbad_in_aper == 0:
			qual[s] = 1
			ngood+=1
		else:
			qual[s] = 0
			nbad+=1


	print("bad: {0}, good: {1} pixels.".format(nbad,ngood))
	flux = np.sum(flux_cube*ix, axis=(1,2))
	##
	## "the variance of a sum of uncorrelated events is the sum of the variances"
	## standard error of the mean = sample standard deviation / sqrt(sample size)
	##
	numpix = np.sum(ix,axis=(1,2))
#	pdb.set_trace()
	noise = np.sqrt(np.sum((noise_cube*ix)**2, axis=(1,2)))/np.sqrt(numpix)
		
	return(flux,noise,qual)


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
	if np.size(dd_pixel) > 1:
		z,y,x=np.mgrid[:flux_cube.shape[0],:flux_cube.shape[1],:flux_cube.shape[2]]
		ix  = np.abs(y-dd_pixel[:,None,None]) < r_px_y
	else:
		z,y,x=np.mgrid[:flux_cube.shape[0],:flux_cube.shape[1],:flux_cube.shape[2]]
		ix  = np.abs(y-dd_pixel) < r_px_y

	##
	##
	flux,noise,qual = qual_interpret_interpolate(flux_cube, noise_cube, qual_cube, ix)
	
	if outfile:		
		## plot (flux vs. wavelength + quality)
		norm=np.median(flux)
		plt.plot(wave,flux/norm)
		plt.plot(wave,noise/norm,'g-')
		plt.plot(wave[qual==0],flux[qual==0]/norm,'rx')
		plt.xlabel("Wavelength")
		plt.ylabel("Flux (normalized)")

		if arm == "UVB":
			plt.ylim([0,3])
		else:
			plt.ylim([0,2])

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
		plt.xlabel("Wavelength")
		plt.ylabel("Flux (normalized)")
		plt.ylim([0,1.5])
		plt.xlim(wregion)
		plt.savefig(outfile+"_"+arm+"_detail.pdf")
		plt.close()
	else:
		plt.show()
		
	return(wave,flux,noise,qual)

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
	dprlist=["SCI","TELL","FLUX"]

	for arm in arms:
		dataset_definition = os.getenv("XDIR")+'/dataset_definition/'+ob_name+'_'+arm+'.txt'
		a=ascii.read(dataset_definition,data_start=0)
				
		for row in a:
			ob=row[0]
			dpid=row[1]
			dpr=dprlist[row.index]
			dpid=dpid.split(".fits")[0]
			dir_cube=os.getenv("XDIRRED")+"/"+dpid.replace(":","_")+"_tpl/"
			f_cube=dir_cube+ob+"_"+dpr+"_IFU_MERGE3D_DATA_OBJ_"+arm+".fits"
			
			if arm=="UVB" and dpr=="TELL":
				continue

			if not os.path.isfile(f_cube):
				raise IOError("Cube file does not exist at", f_cube)

			dir_out=os.getenv("XDIR")+"/spectra/"+ob_name+"/"
			if not os.path.isdir(dir_out):
				os.mkdir(dir_out)
		
			f_out=dir_out+dpr+"_"+arm+"_spec.fits"

			if arm=="NIR" and dpr=="SCI":
				dpid_sci_nir = dpid
				f_cube_sci_nir = f_cube
				f_out_sci_nir = f_out
			if arm=="NIR" and dpr=="TELL":
				dpid_tell_nir = dpid
				f_cube_tell_nir = f_cube
				f_out_tell_nir = f_out
			if arm=="NIR" and dpr=="FLUX":
				dpid_flux_nir = dpid			
				f_cube_flux_nir = f_cube
				f_out_flux_nir = f_out
			if arm=="VIS" and dpr=="SCI":
				dpid_sci_vis = dpid
				f_cube_sci_vis = f_cube
				f_out_sci_vis = f_out
			if arm=="VIS" and dpr=="TELL":
				dpid_tell_vis = dpid
				f_cube_tell_vis = f_cube
				f_out_tell_vis = f_out
			if arm=="VIS" and dpr=="FLUX":
				dpid_flux_vis = dpid
				f_cube_flux_vis = f_cube
				f_out_flux_vis = f_out
			if arm=="UVB" and dpr=="SCI":
				dpid_sci_uvb = dpid
				f_cube_sci_uvb = f_cube
				f_out_sci_uvb = f_out
			if arm=="UVB" and dpr=="FLUX":
				dpid_flux_uvb = dpid
				f_cube_flux_uvb = f_cube
				f_out_flux_uvb = f_out

	##
	## for a few objects automatic extraction fails and we need to set y_NIR manually
	y_SCI_NIR_manually=False
	if dpid_sci_nir=="XSHOO.2015-05-21T07:40:31.524":
		y_SCI_NIR_manually=11.0
	if dpid_sci_nir=="XSHOO.2014-02-21T03:57:08.535": ## NGC3351_1
		y_SCI_NIR_manually=4.0
	if dpid_sci_nir=="XSHOO.2015-05-21T02:55:25.591": ## NGC5128_1
		y_SCI_NIR_manually=9.0

	##
	## do centroiding on near-IR arm and compute atmospheric dispersion correction from that position
	##	
	dd_pixel_SCI_UVB, dd_pixel_SCI_VIS, y_SCI_NIR = dar_position(dpid_sci_nir, DPID_VIS=dpid_sci_vis, DPID_UVB=dpid_sci_uvb, fplot=ob_name+"_SCI",y_NIR=y_SCI_NIR_manually)
	dd_pixel_TELL_UVB, dd_pixel_TELL_VIS, y_TELL_NIR = dar_position(dpid_tell_nir) ## no QC plot here since I normally don't have the UVB arm data for the telluric
	dd_pixel_FLUX_UVB, dd_pixel_FLUX_VIS, y_FLUX_NIR = dar_position(dpid_flux_nir, DPID_VIS=dpid_flux_vis, DPID_UVB=dpid_flux_uvb, fplot=ob_name+"_FLUX")

	flatten_spectrum(f_cube_sci_nir, "NIR", y_SCI_NIR, outfile=f_out_sci_nir)
	flatten_spectrum(f_cube_sci_vis, "VIS", dd_pixel_SCI_VIS, outfile=f_out_sci_vis)
	flatten_spectrum(f_cube_sci_uvb, "UVB", dd_pixel_SCI_UVB, outfile=f_out_sci_uvb)

	## not creating spectrum for UVB arm telluric since it normally is not reduced (since not required)
	flatten_spectrum(f_cube_tell_nir, "NIR", y_TELL_NIR, outfile=f_out_tell_nir)
	flatten_spectrum(f_cube_tell_vis, "VIS", dd_pixel_TELL_VIS, outfile=f_out_tell_vis)

	flatten_spectrum(f_cube_flux_nir, "NIR", y_FLUX_NIR, outfile=f_out_flux_nir)
	flatten_spectrum(f_cube_flux_vis, "VIS", dd_pixel_FLUX_VIS, outfile=f_out_flux_vis)
	flatten_spectrum(f_cube_flux_uvb, "UVB", dd_pixel_FLUX_UVB, outfile=f_out_flux_uvb)

 


################################################################################
#################################### "MAIN" ####################################
################################################################################

#dir="/Volumes/astrodata2/XSHOODATA-Reduced/reflex_end_products/2015-11-12_pipeline_2-6-8/"
#files=ascii.read("/Users/leo/Desktop/files.txt")
#for f in files:
#	ff=dir+f[0]
#	print(ff)
#	w,f,n,q,c = flatten_spectrum(ff)
