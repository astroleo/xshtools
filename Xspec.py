import os
import sqlite3
import numpy as np

from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from astropy.io import ascii
from astropy.table import Table
from astroquery.ned import Ned

from unmask_spectra import unmask_spectra
from postprocess import flatten_ob
from molecfit_wrapper import *

import pdb

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
	def __init__(self, ob_name, object_name=""):
		self.ob_name = ob_name
		self.object_name = object_name

		if self.object_name == "":
			if self.ob_name[:-2] == "ESO137":
				self.object_name = "ESO137-G034"
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
		self.cfg_SL_infile=self.dir_starlight+"/infiles/"+self.ob_name+".in"
		self.dataset_definition = os.getenv("XDIR")+'/dataset_definition/'+self.ob_name+'.txt'
		d=ascii.read(self.dataset_definition)
		self.night=d['night'][0]
		ob_name_list=[d['object'][0],d['telluric'][0],d['flux'][0]]
		self.arms=["NIR","VIS","UVB"]
		self.dprlist=["SCI","TELL","FLUX"]
	
		self.dataset = Table(names=('dprtype', 'ob_name', 'arm', 'dpid'), dtype=(np.dtype((str,10)), np.dtype((str,20)), np.dtype((str,3)), np.dtype((str,29))), meta={'night': self.night})

		conn=sqlite3.connect(os.getenv("OBSDB"))
		c=conn.cursor()
	
		for arm in self.arms:
			f_arm=self.dataset_definition.split('.txt')[0]+'_'+arm+'.txt'
			if os.path.isfile(f_arm):
				d_arm=ascii.read(f_arm,data_start=0,names=["ob","dpid"])
				for ob,dpid,dpr in zip(d_arm["ob"],d_arm["dpid"],self.dprlist):
					self.dataset.add_row([dpr,ob,arm,dpid])
			else:
				with open(f_arm,'w') as f:
					for ob,dpr in zip(ob_name_list,self.dprlist):
						query = "select arcfile from shoot where night=\"" + self.night + \
							"\" and ob_name=\""+ ob + \
							"\" and arm=\"" + arm + \
							"\" and opti2_name=\"IFU\" limit 1;"
						c.execute(query)
						dpid=c.fetchone()
						self.dataset.add_row([dpr,ob,arm,dpid])
						file_string=ob+" "+dpid[0]+"\n"
						f.write(file_string)

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
					for row in self.dataset[np.where(self.dataset['arm']==arm)]:
						dpr = row['dprtype']
						ob = row['ob_name']
						dpid=row['dpid']

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
		flatten_ob(self.ob_name)

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
	def intercalibrate(self,QCplot=False):
		if os.path.isfile(self.f_combined):
			print("intercalibrate: combined FITS file ("+self.f_combined+") already exists. Doing nothing here.")
			return()
		##
		## TODO: check if variables already exist
		self.get_calibrated_data()
		f_intercalibration_QCplot = os.getenv("XDIR") + "/QC/intercalibration/"+self.ob_name+".png"
		##
		## determine VIS correction factor from average flux at end of UVB and beginning of VIS arms
		fu_end = np.median(fu[(qu==1) & (wu > 550) & (wu < 560)])
		fv_start = np.median(fv[(qv==1) & (wv > 560) & (wv < 570)])
		factor_vis = fu_end/fv_start
		##
		## determine NIR correction factor from average flux at end of VIS and beginning of NIR arms
		fv_end = np.median(fv[(qv==1) & (wv > 1000) & (wv < 1020)])
		fn_start = np.median(fn[(qn==1) & (wn > 1020) & (wn < 1040)])
		factor_nir = fv_end/fn_start

		if QCplot and not os.path.isfile(f_intercalibration_QCplot):
			##
			## subplot for UVB - VIS intersection
			plt.subplot(211)
			plt.plot(wv[wv>550],fv[wv>550], color='lightgreen', linestyle="dotted", ms=1)
			plt.plot(wu[wu<560],fu[wu<560], color='lightblue', linestyle="dotted", ms=1)
			plt.xlim([520,600])
			##
			## determine plotting range
			m=np.median(fv[(qv==1) & (wv>560) & (wv<600)])
			s=np.std(fv[(qv==1) & (wv>560) & (wv<600)])
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
			plt.plot(wn[wn>1000],fn[wn>1000],color='peachpuff', linestyle="dotted", ms=1)
			plt.plot(wv[wv<1030],fv[wv<1030], color='lightgreen', linestyle="dotted", ms=1)
			plt.xlim([1000,1050])
			##
			## determine plotting range
			m=np.median(fv[(qv==1) & (wv>1000) & (wv<1050)])
			s=np.std(fv[(qv==1) & (wv>1000) & (wv<1050)])
			plt.ylim([np.min(0,s-3*m),s+5*m])
			##
			## show and print correction factor
			plt.plot([1000,1020],[fv_end,fv_end],'g-',linewidth=2)
			plt.plot([1020,1040],[fn_start,fn_start],'r-',linewidth=2)
			yloc=plt.ylim()[0]+0.1*np.diff(plt.ylim())
			plt.text(1045,yloc,"Correction factor for NIR spectrum (excl. VIS factor): " + str(factor_nir),ha="right")
			plt.suptitle("Flux and telluric calibrated observation for " + ob_name)
			plt.savefig(f_intercalibration_QCplot)
			plt.clf()
		##
		## combine spectra and produce combined outfile
		wave = np.hstack([wu[wu<560],wv[(wv>560) & (wv<1020)],wn[(wn>1020)]])
		flux = np.hstack([fu[wu<560],factor_vis * fv[(wv>560) & (wv<1020)],factor_vis * factor_nir * fn[(wn>1020)]])
		noise = np.hstack([nu[wu<560],factor_vis * nv[(wv>560) & (wv<1020)],factor_vis * factor_nir * nn[(wn>1020)]])
		qual = np.hstack([qu[wu<560],qv[(wv>560) & (qv<1020)],qn[(wn>1020)]])
		
		prihdu=fits.PrimaryHDU(header=hdr_uvb)
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
	
		if os.path.is_file(self.f_starlight_bc03):
			print("write_spec_starlight_bc03: outfile exists in " + self.f_starlight_bc03 + ". Doing nothing.")
			return()
		
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
		## interpolation is not enough; need to also reduce the noise!
		v=interp1d(wave,1/(noise**2))
		ivariance=v(lamgrid)
		inoise=1/np.sqrt(ivariance)
		##
		## reduce noise by 1/sqrt(N) where N is number of bins averaged over
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
		print("Wrote ", outfile)
		
	####################################################################################################
	##
	## METHOD write_SL_infile
	##
	## write STARLIGHT input file and store relevant parameters in data structure of this object for later use, e.g. with plotting functions
	##
	####################################################################################################
	def write_SL_infile(speed="fast"):
		if not os.path.isdir(self.cfg_SL_infile):
			raise ValueError(self.cfg_SL_infile + " does not exist; need to enhance this method to produce it automatically!")

	####################################################################################################
	##
	## METHOD get_SL_params
	##
	## parse STARLIGHT config files and store relevant parameters in data structure of this object for later use, e.g. with plotting functions
	##
	####################################################################################################
	def get_SL_params():
		a=subprocess.run("grep Olsyn_ini " + self.SL_infile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
		Olsyn_ini=np.float(a.stdout.strip())
		a=subprocess.run("grep Olsyn_fin " + self.SL_infile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
		Olsyn_fin=np.float(a.stdout.strip())
		a=subprocess.run("head -n 16 " + self.SL_infile + " | tail -n 1 | awk '{print $2}'", shell=True, universal_newlines=True, stdout=subprocess.PIPE)
		configfile=a.stdout.strip()

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
	def run_starlight():
		cwd=ps.getcwd()
		os.chdir(dir_starlight)
		subprocess.run("./StarlightChains_v04.exe < " + self.SL_infile, shell=True)
		os.chdir(cwd)
		
		Olsyn_ini,Olsyn_fin,llow_norm,lupp_norm = get_SL_params(infile)
		
		subprocess.call(dir_starlight + "/scripts/extract_results.sh",shell=True)

		plotSL(id, pdf=True, Olsyn_ini=Olsyn_ini, Olsyn_fin=Olsyn_fin, llow_norm=llow_norm, lupp_norm=lupp_norm)
		#plot_popvec_hist()
		