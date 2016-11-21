import os
import shutil
import subprocess
import sqlite3
import numpy as np

from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from astropy.io import ascii, fits
from astropy.table import Table
from astroquery.ned import Ned

from unmask_spectra import unmask_spectra
from molecfit_wrapper import *
from AutoVivification import AutoVivification
from dar_extract import dar_position

from SLplotlib import plotSL
from XQC_plot import XQC_plot

import params

import pdb


####################################################################################################
##
## FUNCTION qual_interpret_interpolate
##
## interpret QUAL extension in X-SHOOTER reduced data
##
## INPUT
##     flux_cube, noise_cube
##     qual_cube: cube with bad pixel codes
##     ix: aperture mask
##
####################################################################################################
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
	noise = np.sqrt(np.sum((noise_cube*ix)**2, axis=(1,2)))/np.sqrt(numpix)

	return(flux,noise,qual)







class Xspec():
	"""
	A class to handle X-SHOOTER spectra and provide an interface to kinematic fits (ppxf) and SSP fits (STARLIGHT)
	
	Among other things, this does
	- build file list for a given OB (requires a dataset definition file)
	- extracts spectra from pipeline reduced data with a given aperture and taking into account differential atmospheric refraction
	- flux calibration
	- telluric calibration
	- merging of the different arms ("intercalibration")
	- for each step, "quality control" plots can be produced in order to check what happens
	- final product: telluric-corrected, calibrated, merged spectrum saved to FITS
	- ppxf fit
	- STARLIGHT fit
	"""
	####################################################################################################
	##
	## METHOD __init__
	##
	## build dataset from previously saved files or from queries to $OBSDB if files do not yet exist
	##
	####################################################################################################
	def __init__(self, ob_name, delete_old=False, object_name=""):
		self.ob_name = ob_name
		self.object_name = object_name

		if self.object_name == "":
			if self.ob_name[:-2] == "ESO137":
				self.object_name = "ESO137-G034"
			elif self.ob_name[:-2] == "ESO021":
				self.object_name = "ESO021-G004"
			else:
				self.object_name = ob_name[:-2]

		##
		## TODO: need some error handling here!
		n=Ned.query_object(self.object_name)
		self.z=n['Redshift'].data.data[0]

		##
		## set some directories and files
		self.f_combined = os.getenv("XDIR")+'/combined/'+self.ob_name+'.fits'
		self.dir_starlight = os.getenv("HOME") + "/STARLIGHT"
		self.f_starlight_bc03 = self.dir_starlight + "/spectra/" + self.ob_name + ".txt"
		self.SL_infile=self.dir_starlight+"/"+self.ob_name+".in"
		self.dataset_definition = os.getenv("XDIR")+'/dataset_definition/'+self.ob_name+'.txt'
		d=ascii.read(self.dataset_definition,comment="#")
		self.night=d['night'][0]
		##
		## self.ob_name == self.ob_names["SCI"] -- the reasoning is that self.ob_name is the main reference for a dataset which is referenced by the science OB name
		self.ob_names={"SCI":d['object'][0], "TELL": d['telluric'][0], "FLUX": d['flux'][0]}
		self.arms=["NIR","VIS","UVB"]
		self.dprlist=["SCI","TELL","FLUX"]
		
		##
		## delete old files if requested
		if delete_old:
			caldir = os.getenv("XDIR")+'/calibrated/'+self.ob_name
			molecfit_dir = os.getenv("XDIR")+'/molecfit/'+self.ob_name
			spec_dir = os.getenv("XDIR")+'/spectra/'+self.ob_name
			combined_spec = os.getenv("XDIR")+'/combined/'+self.ob_name+'.fits'
			
			if os.path.isdir(caldir):
				shutil.rmtree(caldir)
				print("Removed directory ", caldir)
			
			if os.path.isdir(molecfit_dir):
				shutil.rmtree(molecfit_dir)
				print("Removed directory ", molecfit_dir)
			
			if os.path.isdir(spec_dir):
				shutil.rmtree(spec_dir)
				print("Removed directory ", spec_dir)
			
			if os.path.isfile(combined_spec):
				os.remove(combined_spec)
				print("Removed file ", combined_spec)
		
		self.dataset = AutoVivification()
		
		self.dir_out=os.getenv("XDIR")+"/spectra/"+self.ob_name+"/"
		if not os.path.isdir(self.dir_out):
			os.mkdir(self.dir_out)

		conn=sqlite3.connect(os.getenv("OBSDB"))
		c=conn.cursor()

		for arm in self.arms:
			f_arm=self.dataset_definition.split('.txt')[0]+'_'+arm+'.txt'
			if not os.path.isfile(f_arm):
				# generate dataset definition file for this arm
				with open(f_arm,'w') as f:
					for dpr in self.dprlist:
						ob=self.ob_names[dpr]
						query = "select arcfile from shoot where night=\"" + self.night + \
							"\" and ob_name=\""+ ob + \
							"\" and arm=\"" + arm + \
							"\" and opti2_name=\"IFU\" limit 1;"
						##
						## TODO -- need some error handling here for the case when the 
						## query does not return any results, e.g. because the
						## dataset definition file is wrong
						##
						c.execute(query)
						dpid=c.fetchone()
						dpid=dpid[0]
						if not dpid.startswith("XSHOO."):
							raise ValueError("Invalid or empty dpid", dpid)
						file_string=ob+" "+dpid+"\n"
						f.write(file_string)
			# parse file and populate self.dataset
			d_arm=ascii.read(f_arm,data_start=0,names=["ob","filename"],comment="#")
			for ob,filename,dpr in zip(d_arm["ob"],d_arm["filename"],self.dprlist):
				dpid = filename.split(".fits")[0]
				self.dataset[dpr][arm]["dpid"] = dpid
				dir_cube = os.getenv("XDIRRED")+"/"+dpid.replace(":","_")+"_tpl/"
				self.dataset[dpr][arm]["dir_cube"] = dir_cube
				f_cube = dir_cube + self.ob_names[dpr] +"_"+dpr+"_IFU_MERGE3D_DATA_OBJ_"+arm+".fits"
				##
				## do not raise error for missing reduced UVB TELL file since it is not required
				if not os.path.isfile(f_cube) and not dpr=="TELL" and not arm=="UVB":
					raise IOError("Cube file does not exist at", f_cube)
				self.dataset[dpr][arm]["f_cube"] = f_cube
				self.dataset[dpr][arm]["f_out"] = self.dir_out + dpr + "_" + arm + "_spec.fits"
				##
				## Caution: f_raw may not be available as it is usually stored on an external disk; 
				##             we check that later when we actually want to read that file
				##
				## determine "night" of this data file from filename
				this_night = subprocess.getoutput("/Users/leo/miditools/f/whichnight_date.sh " + filename[6:25])
				self.dataset[dpr][arm]["f_raw"] = os.getenv("XSHOODATA") + "/" + this_night + "/" + filename				
				self.dir_out + dpr + "_" + arm + "_spec.fits"


	####################################################################################################
	##
	## METHOD flatten_spectrum
	##
	## take data cube, extract 1D spectrum from it, call XQC_plot
	##
	####################################################################################################	
	def flatten_spectrum(self, dpr, arm, dd_pixel):
		f = self.dataset[dpr][arm]["f_cube"]
		if not os.path.isfile(f):
			raise IOError("Reduced data cube file not available: ", f)

		outfile = self.dataset[dpr][arm]["f_out"]
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
		r_px_y = np.ceil(params.aperture_arcsec/(2*scale_y_arcsec))

		##
		## dd_pixel is either an array of same length as wave (check) or a single value (y_NIR from dar_position)
		## here: computing (wavelength-dependent) mask 'ix' from D.A.R. corrected NIR centroid position
		if np.size(dd_pixel) > 1:
			z,y,x=np.mgrid[:flux_cube.shape[0],:flux_cube.shape[1],:flux_cube.shape[2]]
			ix  = np.abs(y-dd_pixel[:,None,None]) < r_px_y
		else:
			z,y,x=np.mgrid[:flux_cube.shape[0],:flux_cube.shape[1],:flux_cube.shape[2]]
			ix  = np.abs(y-dd_pixel) < r_px_y

		flux,noise,qual = qual_interpret_interpolate(flux_cube, noise_cube, qual_cube, ix)
		
		##
		## create QC plot showing raw counts + bad pixel histogram, raw count spectrum + masked slice image
		dataset_raw = self.dataset[dpr][arm]["f_raw"]
		if not os.path.isfile(dataset_raw):
			raise IOError("Raw data not accessible at", dataset_raw)

		###   need to adapt XQC_plot (make part of class?) in order to get access to masked data etc.
		### perhaps it works as such -- i.e. reading data again in XQC_plot...
		XQC_plot(self.dataset[dpr][arm]["dpid"], dataset_raw, f, dpr, self.night, arm, ix, flux_cube, qual_cube)
		XQC_plot(self.dataset[dpr][arm]["dpid"], dataset_raw, f, dpr, self.night, arm, ix, flux_cube, qual_cube, plot_type="full")
		

	
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
		
		return(wave,flux,noise,qual)

	####################################################################################################
	##
	## METHOD flatten_ob
	##
	## prepare all data needed to calibrate one science observation
	##
	####################################################################################################

	def flatten_ob(self):
		##
		## for a few objects automatic extraction fails and we need to set y_NIR manually
		y_SCI_NIR_manually=False
		if self.dataset["SCI"]["NIR"]["dpid"]=="XSHOO.2015-05-21T07:40:31.524":
			y_SCI_NIR_manually=11.0
		if self.dataset["SCI"]["NIR"]["dpid"]=="XSHOO.2014-02-21T03:57:08.535": ## NGC3351_1
			y_SCI_NIR_manually=4.0
		if self.dataset["SCI"]["NIR"]["dpid"]=="XSHOO.2015-05-21T02:55:25.591": ## NGC5128_1
			y_SCI_NIR_manually=9.0

		##
		## do centroiding on near-IR arm and compute atmospheric dispersion correction from that position
		##	
		dd_pixel_SCI_UVB, dd_pixel_SCI_VIS, y_SCI_NIR = dar_position(self.dataset["SCI"]["NIR"]["dpid"], 
			DPID_VIS=self.dataset["SCI"]["VIS"]["dpid"], 
			DPID_UVB=self.dataset["SCI"]["UVB"]["dpid"], 
			fplot=self.ob_name+"_SCI",
			y_NIR=y_SCI_NIR_manually)
		##
		## no QC plot for TELL since I normally don't have the UVB arm data for the telluric
		dd_pixel_TELL_UVB, dd_pixel_TELL_VIS, y_TELL_NIR = dar_position(self.dataset["TELL"]["NIR"]["dpid"])
		dd_pixel_FLUX_UVB, dd_pixel_FLUX_VIS, y_FLUX_NIR = dar_position(self.dataset["FLUX"]["NIR"]["dpid"], 
			DPID_VIS=self.dataset["FLUX"]["VIS"]["dpid"], 
			DPID_UVB=self.dataset["FLUX"]["UVB"]["dpid"], 
			fplot=self.ob_name+"_FLUX")
		
		self.flatten_spectrum("SCI", "NIR", y_SCI_NIR)
		self.flatten_spectrum("SCI", "VIS", dd_pixel_SCI_VIS)
		self.flatten_spectrum("SCI", "UVB", dd_pixel_SCI_UVB)

		## not creating spectrum for UVB arm telluric since it normally is not reduced (since not required)
		self.flatten_spectrum("TELL", "NIR", y_TELL_NIR)
		self.flatten_spectrum("TELL", "VIS", dd_pixel_TELL_VIS)

		self.flatten_spectrum("FLUX", "NIR", y_FLUX_NIR)
		self.flatten_spectrum("FLUX", "VIS", dd_pixel_FLUX_VIS)
		self.flatten_spectrum("FLUX", "UVB", dd_pixel_FLUX_UVB)

	####################################################################################################
	##
	## METHOD build_dataset_QC_page
	##
	## build HTML QC pages by dataset
	##
	####################################################################################################
	def build_dataset_QC_page(self):
		html_QC = os.getenv("XDIR")+'/QC/datasets/'+self.ob_name+'.html'
		if os.path.isfile(html_QC):
			print(html_QC + ' exists.')
		else:
			with open(html_QC,'w') as html:
				html.write("<html>\n")
				html.write("<header>\n")
				html.write("<title>"+self.ob_name+"</title>\n")
				html.write("</header>\n")
				html.write("<body>\n")
				
				for arm in self.arms:
					html.write("<h1>"+arm+"</h1>")
					for dpr in self.dprlist:
						ob = self.ob_names[dpr]
						dpid = self.dataset[dpr][arm]["dpid"]

						html.write("<h2>" + dpr + ": " + ob + "</h2>\n")
						img="../plot/" + self.night + "/" + dpid.replace(":","_") + "_" + dpr + "_" + arm + ".png"
						html.write("<img src=\"" + img + "\">\n")
						html.write("<br>\n")
						html.write("<hr>\n")

				html.write("</body>\n")
				html.write("</html>\n")

	####################################################################################################
	##
	## METHOD reduce_data
	##
	## wrapper function to organize the data reduction, flux calibration, telluric correction, ...
	##
	####################################################################################################
	def reduce_data(self):
		self.flatten_ob()

		run_molecfit(self.ob_name,"NIR")
		molecfit_QC(self.ob_name,"NIR")
		run_molecfit(self.ob_name,"VIS")
		molecfit_QC(self.ob_name,"VIS")

		flux_calibrate(self.ob_name,"NIR")
		flux_calibrate(self.ob_name,"VIS")
		flux_calibrate(self.ob_name,"UVB")

	
	####################################################################################################
	##
	## METHOD get_calibrated_data
	##
	## does what it says (and also checks whether data exists)
	##
	####################################################################################################
	def get_calibrated_data(self,reduce=False):
		dir=os.getenv("XDIR")+"/calibrated/"+self.ob_name+"/"
		self.uvbfile=dir+"SCI_UVB_calibrated.fits"
		self.visfile=dir+"SCI_VIS_calibrated.fits"
		self.nirfile=dir+"SCI_NIR_calibrated.fits"
		
		if not os.path.isfile(self.uvbfile):
			raise IOError("calibrated UVB spectrum (" + self.uvbfile + ") does not exist.")
			reduce=True
		if not os.path.isfile(self.visfile):
			raise IOError("calibrated VIS spectrum (" + self.visfile + ") does not exist.")
			reduce=True
		if not os.path.isfile(self.nirfile):
			raise IOError("calibrated NIR spectrum (" + self.nirfile + ") does not exist.")
			reduce=True

		if reduce==True:
			self.reduce_data()
	
		self.hdu_uvb=fits.open(self.uvbfile)
		self.hdu_vis=fits.open(self.visfile)
		self.hdu_nir=fits.open(self.nirfile)

		self.tu=self.hdu_uvb[1].data
		self.hdr_uvb = self.hdu_uvb[0].header
		self.wu=self.tu['WAVE']
		self.fu=self.tu['FLUX']
		self.nu=self.tu['NOISE']
		self.qu=self.tu['QUAL']

		self.tv=self.hdu_vis[1].data
		self.hdr_vis = self.hdu_vis[0].header
		self.wv=self.tv['WAVE']
		self.fv=self.tv['FLUX']
		self.nv=self.tv['NOISE']
		self.qv=self.tv['QUAL']

		self.tn=self.hdu_nir[1].data
		self.hdr_nir = self.hdu_nir[0].header
		self.wn=self.tn['WAVE']
		self.fn=self.tn['FLUX']
		self.nn=self.tn['NOISE']
		self.qn=self.tn['QUAL']

	####################################################################################################
	##
	## METHOD intercalibrate
	##
	## merge UVB-VIS-NIR arms, find flux correction factors, produce QC plot if specified, produce final merged spectrum if it does not exist
	##
	####################################################################################################
	def intercalibrate(self,QCplot=True):
		if os.path.isfile(self.f_combined):
			print("intercalibrate: overwriting existing FITS file ("+self.f_combined+").")
		##
		## TODO: check if variables already exist
		self.get_calibrated_data()
		f_intercalibration_QCplot = os.getenv("XDIR") + "/QC/intercalibration/"+self.ob_name+".png"
		##
		## determine VIS correction factor from average flux at end of UVB and beginning of VIS arms
		fu_end = np.median(self.fu[(self.qu==1) & (self.wu > 550) & (self.wu < 560)])
		fv_start = np.median(self.fv[(self.qv==1) & (self.wv > 560) & (self.wv < 570)])
		factor_vis = fu_end/fv_start
		##
		## determine NIR correction factor from average flux at end of VIS and beginning of NIR arms
		fv_end = np.median(self.fv[(self.qv==1) & (self.wv > 1000) & (self.wv < 1020)])
		fn_start = np.median(self.fn[(self.qn==1) & (self.wn > 1020) & (self.wn < 1040)])
		factor_nir = fv_end/fn_start

		if QCplot and not os.path.isfile(f_intercalibration_QCplot):
			##
			## subplot for UVB - VIS intersection
			plt.subplot(211)
			plt.plot(self.wv[self.wv>550],self.fv[self.wv>550], color='lightgreen', linestyle="dotted", ms=1)
			plt.plot(self.wu[self.wu<560],self.fu[self.wu<560], color='lightblue', linestyle="dotted", ms=1)
			plt.xlim([520,600])
			##
			## determine plotting range
			m=np.median(self.fv[(self.qv==1) & (self.wv>560) & (self.wv<600)])
			s=np.std(self.fv[(self.qv==1) & (self.wv>560) & (self.wv<600)])
			plt.ylim([np.min(0,s-3*m),s+3*m])
			##
			## show and print correction factor
			plt.plot([550,560],[fu_end,fu_end],'b-',linewidth=2)
			plt.plot([560,570],[fv_start,fv_start],'g-',linewidth=2)
			yloc=plt.ylim()[0]+0.1*np.diff(plt.ylim())
			plt.text(595,yloc,"Correction factor for VIS spectrum: " + str(factor_vis),ha="right")
			##
			## subplot for VIS - NIR intersection
			plt.subplot(212)
			plt.plot(self.wn[self.wn>1000],self.fn[self.wn>1000],color='peachpuff', linestyle="dotted", ms=1)
			plt.plot(self.wv[self.wv<1030],self.fv[self.wv<1030], color='lightgreen', linestyle="dotted", ms=1)
			plt.xlim([1000,1050])
			##
			## determine plotting range
			m=np.median(self.fv[(self.qv==1) & (self.wv>1000) & (self.wv<1050)])
			s=np.std(self.fv[(self.qv==1) & (self.wv>1000) & (self.wv<1050)])
			plt.ylim([np.min(0,s-3*m),s+5*m])
			##
			## show and print correction factor
			plt.plot([1000,1020],[fv_end,fv_end],'g-',linewidth=2)
			plt.plot([1020,1040],[fn_start,fn_start],'r-',linewidth=2)
			yloc=plt.ylim()[0]+0.1*np.diff(plt.ylim())
			plt.text(1045,yloc,"Correction factor for NIR spectrum (excl. VIS factor): " + str(factor_nir),ha="right")
			plt.suptitle("Flux and telluric calibrated observation for " + self.ob_name)
			plt.savefig(f_intercalibration_QCplot)
			plt.clf()
		##
		## combine spectra and produce combined outfile
		wave = np.hstack([self.wu[self.wu<560],self.wv[(self.wv>560) & (self.wv<1020)],self.wn[(self.wn>1020)]])
		flux = np.hstack([self.fu[self.wu<560],factor_vis * self.fv[(self.wv>560) & (self.wv<1020)],factor_vis * factor_nir * self.fn[(self.wn>1020)]])
		noise = np.hstack([self.nu[self.wu<560],factor_vis * self.nv[(self.wv>560) & (self.wv<1020)],factor_vis * factor_nir * self.nn[(self.wn>1020)]])
		qual = np.hstack([self.qu[self.wu<560],self.qv[(self.wv>560) & (self.qv<1020)],self.qn[(self.wn>1020)]])
		
		prihdu=fits.PrimaryHDU(header=self.hdr_uvb)
		##
		## TODO: need to remove UVB specific cards!
		col1=fits.Column(name="WAVE", format='1E', array=wave)
		col2=fits.Column(name="FLUX", format='1E', array=flux)
		col3=fits.Column(name="NOISE", format='1E', array=noise)
		col4=fits.Column(name="QUAL", format='1J', array=qual)
		cols=fits.ColDefs([col1,col2,col3,col4])
		tbhdu=fits.BinTableHDU.from_columns(cols)
		tbhdulist = fits.HDUList([prihdu, tbhdu])
		tbhdulist.writeto(self.f_combined,clobber=True)
	
	####################################################################################################
	##
	## METHOD write_spec_starlight_bc03
	##
	## smoothes spectrum to BC03 resolution and writes it out into STARLIGHT format
	##
	## arguments:
	##    wrange         restframe wavelength range to use (must obviously be smaller 
	##                     than de-redshifted wavelength range of input spectrum)
	##                     wrange is given in Angstrom, inclusive range, i.e. start *and* end points are included
	##
	####################################################################################################
	def write_spec_starlight_bc03(self,wrange=[3800,10000],QCplot=False):
		##
		## ensure that relevant data files are there
		self.intercalibrate()
	
		if os.path.isfile(self.f_starlight_bc03):
			print("write_spec_starlight_bc03: Overwriting existing file in " + self.f_starlight_bc03)
		
		wave,flux,noise = unmask_spectra(self.f_combined)
		wave*=10 ## in Angstrom
		##
		## de-redshift spectrum
		wave/=(1+self.z)
		##
		## check that input wavelength range is larger than output wavelength range
		if (wave[0] > wrange[0]) or (wave[-1] < wrange[-1]):
			raise ValueError("input wavelength range must be larger than output wavelength range")
		##
		## normalize flux and noise
		m=np.median(flux)
		flux/=m
		noise/=m
		##
		## define output wavelength vector
		wsampling=1 ## wavelength sampling in Angstrom
		lamgrid=wrange[0] + wsampling * np.arange(1+np.floor((wrange[1]-wrange[0])/wsampling))
		##
		## convolve to resolution of spectral base, i.e. BC03 res = 3 \AA
		##
		## XSHOO resolution is 12600 (VIS) independent of lambda, i.e. choose kernel 
		##    appropriate for central wavelength
		xres=12600
		FWHM_gal=np.average(wrange)/xres
		FWHM_tem=3
		FWHM_dif = np.sqrt(FWHM_tem**2 - FWHM_gal**2)
		sigma = FWHM_dif/2.355/(wave[2]-wave[1]) # Sigma difference in pixels
		flux=gaussian_filter1d(flux,sigma)
		##
		## interpolation is flux-conserving when accounting for different spectral pixel size
		##    dl is 0.2 in UVB and VIS, 0.6 in NIR
		##
		f=interp1d(wave,flux)
		iflux=f(lamgrid)
		##
		## interpolate inverse variance (since this is the additive quantity)
		##
		v=interp1d(wave,1/(noise**2))
		ivariance=v(lamgrid)
		inoise=1/np.sqrt(ivariance)
		##
		## reduce noise by 1/sqrt(N) where N is number of bins averaged
		N = wsampling/np.diff(wave)
		N2=np.hstack([N[0],N]) ## np.diff reduces size of array by 1
		Ni = interp1d(wave,N2)
		Ninterp = Ni(lamgrid)
		fac=np.sqrt(Ninterp)
		inoise/=fac
		## masks for transition regions of XSHOOTER spectral arms
		imask=np.zeros(iflux.shape)
	#	add_mask_telluric(lamgrid, imask, [5530,5600], z)
	#	add_mask_telluric(lamgrid, imask, [10100,10200], z)

		## Quality Check plot
		##
		if QCplot:
			plt.plot(wave,flux,label="intrinsic resolution")
			plt.plot(lamgrid,iflux,label="downsampled resolution")
			##
			## pick region around CaT if within region
			CaT1=8400
			CaT2=8800
	
			if (8000 > wrange[0]) & (9000 < wrange[1]):
				plt.xlim([8000,9000])
				plt.ylim([0.8,1.2])
			else:
				plt.xlim(wrange)
				plt.ylim([0,2])

			plt.legend(loc=2)
			outfile_qc=self.f_starlight_bc03.split('.txt')[0]+'_interpolate_QC.pdf'
			plt.savefig(outfile_qc)
			print("Saved QC plot to ", outfile_qc)
	
		np.savetxt(self.f_starlight_bc03,np.transpose((lamgrid,iflux,inoise,imask)),fmt=('%5.0f.   %8.3f   %8.3f    %3.0f'))
		print("Wrote ", self.f_starlight_bc03)
		
	####################################################################################################
	##
	## METHOD write_SL_infile
	##
	## write STARLIGHT input file and store relevant parameters in data structure of this object for later use, e.g. with plotting functions
	##
	####################################################################################################
	def write_SL_infile(self,speed="fast"):
		if not os.path.isdir(self.SL_infile):
			f_template=os.getenv("HOME")+"/STARLIGHT/STARLIGHT.in.template"
			shutil.copy(f_template,self.SL_infile)
			
			with open(self.SL_infile,'a') as f:
				f.write(self.ob_name+".txt   StCv04.C11.config   Base.BC03.N   "+self.ob_name+"   CAL   31   206.0   " + self.ob_name+".out\n")
			
			print("Wrote " + self.SL_infile)
	
	####################################################################################################
	##
	## METHOD get_SL_params
	##
	## parse STARLIGHT config files and store relevant parameters in data structure of this object for later use, e.g. with plotting functions
	##
	####################################################################################################
	def get_SL_params(self):
		a=subprocess.run("grep Olsyn_ini " + self.SL_infile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
		Olsyn_ini=np.float(a.stdout.strip())
		a=subprocess.run("grep Olsyn_fin " + self.SL_infile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
		Olsyn_fin=np.float(a.stdout.strip())
		a=subprocess.run("head -n 16 " + self.SL_infile + " | tail -n 1 | awk '{print $2}'", shell=True, universal_newlines=True, stdout=subprocess.PIPE)
		configfile=self.dir_starlight+"/"+a.stdout.strip()

		a=subprocess.run("grep llow_norm " + configfile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
		llow_norm=np.float(a.stdout.strip())
		a=subprocess.run("grep lupp_norm " + configfile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
		lupp_norm=np.float(a.stdout.strip())
		return(Olsyn_ini,Olsyn_fin,llow_norm,lupp_norm)
	
	####################################################################################################
	##
	## METHOD run_starlight
	##
	## run STARLIGHT with pre-configured configuration
	##
	####################################################################################################	
	def run_starlight(self):
		cwd=os.getcwd()
		os.chdir(self.dir_starlight)
		subprocess.run("./StarlightChains_v04.exe < " + self.SL_infile, shell=True)
		os.chdir(cwd)
		
		##
		## needs some error handling here to check if STARLIGHT has produced all necessary files
		
		Olsyn_ini,Olsyn_fin,llow_norm,lupp_norm = self.get_SL_params()
		
		subprocess.call(self.dir_starlight + "/scripts/extract_results.sh",shell=True)

		plotSL(self.ob_name, self.SL_infile, pdf=True, Olsyn_ini=Olsyn_ini, Olsyn_fin=Olsyn_fin, llow_norm=llow_norm, lupp_norm=lupp_norm)
		#plot_popvec_hist()
		