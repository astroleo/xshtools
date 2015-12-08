import os
import shutil
from astropy.io import ascii, fits
import numpy as np
from matplotlib import pyplot as plt
import subprocess
from scipy.interpolate import interp1d

##
## generate molecfit config files for correcting an OB as specified in a
##    dataset_definition file and run molecfit
##
def run_molecfit(ob_name):
	arm="VIS"
	xdir=os.getenv('PROJECTS')+'/LP-BAT/XSHOOTER/'
	mtemplate=xdir+'data/molecfit/molecfit_VIS_template.par'
	mdir=xdir+'data/molecfit/'+ob_name
	if not os.path.isdir(mdir):
		os.mkdir(mdir)
		
	dprlist=["SCI","TELL","FLUX"]
	
	spec_sci=xdir+"data/spectra/"+ob_name+"/SCI_VIS_spec.fits"
	if not os.path.isfile(spec_sci):
		print("File missing:",spec_sci)
		return()

	spec_tell=xdir+"data/spectra/"+ob_name+"/TELL_VIS_spec.fits"
	if not os.path.isfile(spec_tell):
		print("File missing:",spec_tell)
		return()

	spec_flux=xdir+"data/spectra/"+ob_name+"/FLUX_VIS_spec.fits"
	if not os.path.isfile(spec_flux):
		print("File missing:",spec_flux)
		return()

	##
	## generate config files for telluric correction of science object
	##    i.e. telluric correction of telluric + application on science object
	##
	mlist=mdir+"/"+arm+"_TELL.list"
	if os.path.isfile(mlist):
		print("List file exists:",mlist)
		return()
	with open(mlist,"a") as f:
		f.write(spec_tell+"\n")
		f.write(spec_sci+"\n")

	mconfig=mdir+'/'+arm+'_TELL.par'
	if os.path.isfile(mconfig):
		print("Config file exists:",mconfig)
		return()
	shutil.copyfile(mtemplate,mconfig)
	with open(mconfig,"a") as f:
		f.write("filename: "+spec_tell+"\n")
		f.write("listname: "+mlist+"\n")
		f.write("output_dir: "+mdir+"\n")
		f.write("output_name: spec_tell\n")
		f.write("end\n")
	
	subprocess.call("/Users/leo/molecfit/bin/molecfit " + mconfig, 
		shell=True)
	subprocess.call("/Users/leo/molecfit/bin/calctrans " + mconfig, 
		shell=True)
	subprocess.call("/Users/leo/molecfit/bin/corrfilelist " + mconfig, 
		shell=True)
	
	##
	## generate config files for telluric correction of flux standard
	##
	mlist=mdir+"/"+arm+"_FLUX.list"
	if os.path.isfile(mlist):
		print("List file exists:",mlist)
		return()
	with open(mlist,"a") as f:
		f.write(spec_flux+"\n")

	mconfig=mdir+'/'+arm+'_FLUX.par'
	if os.path.isfile(mconfig):
		print("Config file exists:",mconfig)
		return()
	shutil.copyfile(mtemplate,mconfig)
	with open(mconfig,"a") as f:
		f.write("filename: "+spec_flux+"\n")
		f.write("listname: "+mlist+"\n")
		f.write("output_dir: "+mdir+"\n")
		f.write("output_name: spec_flux\n")
		f.write("end\n")
	
	subprocess.call("/Users/leo/molecfit/bin/molecfit " + mconfig, 
		shell=True)
	subprocess.call("/Users/leo/molecfit/bin/calctrans " + mconfig, 
		shell=True)
	subprocess.call("/Users/leo/molecfit/bin/corrfilelist " + mconfig, 
		shell=True)

##
## plot TAC correction for science spectrum
##
def molecfit_QC_plot(spec_TAC):
	if not os.path.isfile(spec_TAC):
		print("File missing:",spec_TAC)
		return()

	f_plot=spec_TAC.split(".")[0]+'.png'
	hdu=fits.open(spec_TAC)
	t=hdu[1].data
	w=t['WAVE']
	f=t['FLUX']
	q=t['QUAL']
	tac_f=t['tacflux']
	tac_df=t['tacdflux']
	##
	## tac_q is different than q (perhaps where the molecfit was bad?)
	## Here I combine both quality indicators, i.e. require a spectral point to be good both in the original spectrum and in the correction spectrum
	tac_q=t['tacqual']
	q_combined=np.all([q,tac_q],axis=0)
	plt.plot(w,f,'magenta')
	plt.plot(w,tac_f,'k')
	plt.plot(w,tac_df,'grey')
	plt.plot(w[q_combined==0],tac_f[q_combined==0],'rx')
	m=np.median(tac_f)
	s=np.std(tac_f)
	plt.ylim([0,m+3*s])
	plt.xlabel("Wavelength [nm]")
	plt.ylabel("ADU")
	tstring=spec_TAC.split("/")[-1].split("_")[0]+" observation for OB  "+spec_TAC.split("/")[-2]
	plt.title(tstring)
	plt.savefig(f_plot)
	
	plt.xlim([800,900])
	plt.savefig(f_plot.split(".")[0]+"_detail_800.png")
	plt.xlim([900,1000])
	plt.savefig(f_plot.split(".")[0]+"_detail_900.png")
	plt.clf()

##
## plot quality of telluric fit for telluric / flux standard star
##
def molecfit_QC_fit_plot(spec_model):
	if not os.path.isfile(spec_model):
		print("File missing:",spec_model)
		return()
	
	f_plot1=spec_model.split(".")[0]+'_1.png'
	f_plot2=spec_model.split(".")[0]+'_2.png'
	hdu=fits.open(spec_model)
	t=hdu[1].data
	w=t['lambda']
	f_obs=t['flux']
	f_model=t['mflux']
	
	plt.plot(w,(f_obs-f_model)/f_obs)
	plt.xlabel("Wavelength [nm]")
	plt.ylabel("(F_obs - F_model)/F_obs")
	plt.title(spec_model.split("/")[-1].split("_")[1] + " observation for OB " +   spec_model.split("/")[-2])
	plt.ylim([-0.1,0.1])
	plt.xlim([0.7615,0.7705])
	plt.savefig(f_plot1)

	plt.xlim([0.939,0.952])
	plt.savefig(f_plot2)
	plt.clf()
	
##
## generate some analysis plots to evaluate quality of the molecfit correction
##
def molecfit_QC(ob_name):
	xdir=os.getenv('PROJECTS')+'/LP-BAT/XSHOOTER/'
	mdir=xdir+'data/molecfit/'+ob_name

	spec_sci_TAC = mdir+"/SCI_VIS_spec_TAC.fits"
	molecfit_QC_plot(spec_sci_TAC)

	spec_tell_TAC = mdir+"/TELL_VIS_spec_TAC.fits"
	molecfit_QC_plot(spec_tell_TAC)

	spec_flux_TAC = mdir+"/FLUX_VIS_spec_TAC.fits"
	molecfit_QC_plot(spec_flux_TAC)
	
	spec_tell_model = mdir+"/spec_tell_fit.fits"
	molecfit_QC_fit_plot(spec_tell_model)

	spec_flux_model = mdir+"/spec_flux_fit.fits"
	molecfit_QC_fit_plot(spec_flux_model)

##
## helper function to get appropriate model spectrum for given flux standard
##    and interpolate to wavelength grid of observed spectrum
##
def interpol_fluxstd_modelspec(spec_flux_TAC,wave):
	xdir=os.getenv('PROJECTS')+'/LP-BAT/XSHOOTER/'

	hdu=fits.open(spec_flux_TAC)
	flux_name=hdu[0].header['HIERARCH ESO OBS NAME'].split("_")[0]
	##
	## WARNING: arm is hardcoded here!
	flux_model_spec=xdir+"data/specphot_VIS/"+flux_name+"_VIS.txt"
	d=np.genfromtxt(flux_model_spec)
	w=d[:,0]
	f=d[:,1]
	interpol_function=interp1d(w,f)
	f_model_interpol=interpol_function(wave)
	return(f_model_interpol,flux_name)

##
## flux-calibrate science spectrum, generate combined mask, FITS files, QC plots
##
def flux_calibrate(ob_name):
	xdir=os.getenv('PROJECTS')+'/LP-BAT/XSHOOTER/'
	mdir=xdir+'data/molecfit/'+ob_name
	caldir=xdir+'data/calibrated/'+ob_name
	if not os.path.isdir(caldir):
		os.mkdir(caldir)

	spec_sci_TAC = mdir+"/SCI_VIS_spec_TAC.fits"
	if not os.path.isfile(spec_sci_TAC):
		print("File missing:",spec_sci_TAC)
#		return()

	spec_flux_TAC = mdir+"/FLUX_VIS_spec_TAC.fits"
	if not os.path.isfile(spec_flux_TAC):
		print("File missing:",spec_flux_TAC)
#		return()
	
	hdu_sci=fits.open(spec_sci_TAC)
	hdr=hdu_sci[0].header
	t_sci=hdu_sci[1].data
	w=t_sci['WAVE']
	q_sci=t_sci['QUAL']
	tac_f_sci=t_sci['tacflux']
	tac_df_sci=t_sci['tacdflux']
	tac_q_sci=t_sci['tacqual']
	q_combined_sci=np.all([q_sci,tac_q_sci],axis=0)
	
	hdu_flux=fits.open(spec_flux_TAC)
	hdu_flux[0].header
	t_flux=hdu_flux[1].data
	q_flux=t_flux['QUAL']
	tac_f_flux=t_flux['tacflux']
	tac_q_flux=t_flux['tacqual']
	q_combined=np.all([q_combined_sci,q_flux,tac_q_flux],axis=0)
	
	f_model_interpol, flux_name = interpol_fluxstd_modelspec(spec_flux_TAC,w)

	plt.plot(w,tac_f_flux/np.median(tac_f_flux),'magenta',label="observed (TAC)")
	plt.plot(w,f_model_interpol/np.median(f_model_interpol),'k',label="model")
	plt.title("Flux STD observation of " + flux_name + "(" + hdu_flux[0].header["DATE-OBS"] + ")")
	plt.xlabel("Wavelength [nm]")
	plt.ylabel("flux/median(flux)")
	plt.savefig(mdir+"/FLUX_VIS_obs_model.png")
	plt.clf()

	flux_corr_factor = f_model_interpol/tac_f_flux
	f_sci_tac_flux = flux_corr_factor * tac_f_sci
	df_sci_tac_flux = flux_corr_factor * tac_df_sci
	
	plt.plot(w,f_sci_tac_flux)
	plt.plot(w,df_sci_tac_flux,'grey')
	plt.plot(w[q_combined==0],f_sci_tac_flux[q_combined==0],'rx')
	plt.xlabel("Wavelength [nm]")
	plt.ylabel(r"$F_{\lambda} [erg/(s \cdot cm^2 \AA)]$")
	plt.title("Flux and telluric calibrated spectrum of " + ob_name,fontsize=12)
	m=np.median(f_sci_tac_flux)
	s=np.std(f_sci_tac_flux)
	plt.ylim([0,m+s])
	plt.savefig(caldir+"/SCI_VIS_calibrated.png")
	
	plt.xlim([600,700])
	plt.savefig(caldir+"/SCI_VIS_calibrated_detail600.png")
	plt.xlim([800,900])
	plt.savefig(caldir+"/SCI_VIS_calibrated_detail800.png")
	
	##
	## store all relevant data in FITS file
	##
	outfile=caldir+"/SCI_VIS_calibrated.fits"
	prihdu=fits.PrimaryHDU(header=hdr)
	col1=fits.Column(name="WAVE", format='1E', array=w)
	col2=fits.Column(name="FLUX", format='1E', array=f_sci_tac_flux)
	col3=fits.Column(name="NOISE", format='1E', array=df_sci_tac_flux)
	col4=fits.Column(name="QUAL", format='1J', array=q_combined)
	cols=fits.ColDefs([col1,col2,col3,col4])
	tbhdu=fits.BinTableHDU.from_columns(cols)
	tbhdulist = fits.HDUList([prihdu, tbhdu])
	tbhdulist.writeto(outfile,clobber=True)

#ob_name="ESO208-G021_1"
#run_molecfit(ob_name)
#molecfit_QC(ob_name)
#flux_calibrate(ob_name)

#ob_name="NGC1079_1"
#run_molecfit(ob_name)
#molecfit_QC(ob_name)
#flux_calibrate(ob_name)

ob_name="NGC6814_1"
run_molecfit(ob_name)
molecfit_QC(ob_name)
flux_calibrate(ob_name)
