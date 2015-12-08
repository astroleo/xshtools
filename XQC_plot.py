# QC plot for each obs / raw data
#- dataset name
#- hist of raw data counts
#- plot of spectrum
#- total number of bad pixels

##
## produce a Quality Check plot for an X-SHOOTER reduced data set
##
import matplotlib
matplotlib.use('agg')

import numpy as np
from astropy.io import ascii,fits
from matplotlib import pyplot as plt
import pdb
from astropy.table import Table
import os
import subprocess
import time
from postprocess import aperture_mask


#dir_reduced="/Users/leo/Projekte/LP-BAT/XSHOOTER/data/reflex_end_products/2015-11-12_pipeline_2-6-8/"
#dir_raw="/Volumes/astrodata2/XSHOODATA/"

#dataset_id="XSHOO.2015-05-13T02_10_54.281"
#dataset_id_eso="XSHOO.2015-05-13T02:10:54.281"
#dataset_red = dir_reduced + dataset_id +  "_tpl/Hip072154_TELL_IFU_MERGE3D_DATA_OBJ_VIS.fits"
#dataset_raw = dir_raw + "2015-05-12/" + dataset_id_eso + ".fits"

#dataset_id_eso="XSHOO.2015-05-13T10:51:59.198"
#dataset_raw="/Volumes/astrodata2/XSHOODATA/2015-05-12/XSHOO.2015-05-13T10:51:59.198.fits"
#dataset_red="/Users/leo/Projekte/LP-BAT/XSHOOTER/data/reflex_end_products/2015-11-12_pipeline_2-6-8/XSHOO.2015-05-13T10_51_59.198_tpl/Calibration_TELL_IFU_MERGE3D_DATA_OBJ_VIS.fits"

#dataset_id_eso="XSHOO.2015-05-13T09:56:59.151"
#dataset_raw="/Volumes/astrodata2/XSHOODATA/2015-05-12/XSHOO.2015-05-13T09:56:59.151.fits"
#dataset_red="/Users/leo/Projekte/LP-BAT/XSHOOTER/data/reflex_end_products/2015-11-12_pipeline_2-6-8/XSHOO.2015-05-13T09_56_59.151_tpl/LTT7987_onoff_IFU_TELL_IFU_MERGE3D_DATA_OBJ_VIS.fits"


def XQC_plot(dataset_id_eso, dataset_raw, dataset_red, dprtype, night, arm):
	## check if required data exist
	if not os.path.isfile(dataset_raw):
		print("RAW file missing: ", dataset_raw)
		return()
	if not os.path.isfile(dataset_red):
		print("REDUCED file missing: ", dataset_red)
		return()

	## check if plot already exists
	dir_out="QC/plot/"+night
	if not os.path.isdir(dir_out):
		os.mkdir(dir_out)
	f_out=dir_out+"/"+dataset_id_eso.replace(":","_")+"_"+dprtype+"_"+arm+".png"
	if os.path.isfile(f_out):
		print("QC PLOT already exists: ", f_out)
		return()

	## data preparation
	hdu_raw=fits.open(dataset_raw)
	d=hdu_raw[0].data
	ob_name=hdu_raw[0].header["HIERARCH ESO OBS NAME"]
	am=hdu_raw[0].header["HIERARCH ESO TEL AIRM START"]
	fwhm=hdu_raw[0].header["HIERARCH ESO TEL AMBI FWHM START "]

	ix,ix2,wave,flux_cube,noise_cube,qual_cube = aperture_mask(dataset_red)
	
	##
	## mask qual and flux cubes
	m_qual_cube=ix*qual_cube
	m_flux_cube=ix*flux_cube
	
	qual_spec=np.any(m_qual_cube!=0,axis=(1,2))

	## spec extraction
	img_spec=np.sum(m_flux_cube,axis=(1,2))
	x=np.arange(len(img_spec))
	
	x1=1000
	x2=24000
	
	x1_detail=15000
	x2_detail=16000
	
	if arm=="UVB":
		x2=13500
		x1_detail=10000
		x2_detail=11000
	
	### BEGIN OF PLOT ###

	# hist of raw data counts
	plt.subplot2grid((3,6),(0,0),colspan=3)
	plt.hist(d.ravel(),bins=50)
	plt.yscale("log")
	plt.xlim([0,2**16-1])
	plt.xticks(20000*np.arange(4))
	plt.title("Histogram of values in RAW data")

	# hist of bad pixel values
	plt.subplot2grid((3,6),(0,3),colspan=3)
	plt.hist(m_qual_cube.ravel(),bins=128,range=[0,127])
	hist_str=str("Total number of bad slices: {0}").format(np.sum(qual_spec))
	plt.yscale("log")
	plt.ylim([1,10**7])
	plt.text(127,10**6,hist_str,ha="right",size=8)
	plt.title("Histogram of bad pixel values")

	# plot of entire spectrum
	plt.subplot2grid((3,6),(1,0),colspan=5)
	plt.plot(x[x1:x2],img_spec[x1:x2])
	m=np.median(img_spec[x1:x2])
	s=np.std(img_spec[x1:x2])
	plt.ylim([m-3*s,m+3*s])
	plt.xlabel("Slice")
	plt.ylabel("ADU")
	title_str=str("{0} (seeing: {1}\", airmass: {2})").format(ob_name,fwhm,am)
	plt.title(title_str)

	# plot of image
	plt.subplot2grid((3,6),(1,5))
	img_collapsed=np.sum(m_flux_cube,axis=0)
	plt.imshow(img_collapsed,origin="lower")
	plt.axis("off")



	# plot of spectrum and bad pixels (integrated over entire cube)
	plt.subplot2grid((3,6),(2,0),colspan=5)
	plt.plot(x,img_spec)
	plt.plot(x[qual_spec],img_spec[qual_spec],'rx')
	plt.xlim([x1_detail,x2_detail])
	m=np.median(img_spec[x1_detail:x2_detail])
	s=np.std(img_spec[x1_detail:x2_detail])
	plt.ylim([m-3*s,m+3*s])
	plt.xlabel("Slice")
	plt.ylabel("ADU")


	# image in a couple different wavelength ranges
	plt.subplot2grid((3,6),(2,5))
	img1=np.sum(m_flux_cube[x1:x2],axis=(0))
	plt.imshow(img1,origin="lower")
	plt.axis("off")

	plt.tight_layout()

	#plt.suptitle(dataset_id_eso)
	#plt.show()
		
	plt.savefig(f_out)
	print("Saved ", f_out)
	plt.close()


def QC_night(night):
	print()
	print()
	print("Now entering night: ", night)
	print()
	print()

	dir_out="QC/plot/"+night
	if not os.path.isdir(dir_out):
		os.mkdir(dir_out)
	
	f_obslist="QC/plot/"+night+"/obslist.txt"
	if not os.path.isfile(f_obslist):
		print("obslist is missing for this night. Constructing...")
		subprocess.Popen("bash $XSHOOTOOLS/ob_paths_night.sh " + night + " > " + f_obslist ,shell=True)
		time.sleep(10)
	
#	t_sleep=0
#	while not os.path.isfile(f_obslist):
#		print("Waiting for ob_paths_night.sh to finish...")
#		time.sleep(1)
#		t_sleep+=1
#		if t_sleep > 10:
#			print("I am tired. And something went wrong.")
	
	obslist=Table.read(f_obslist, format="ascii", data_start=0)

	for obs in obslist:
		night=obs[0]
		arm=obs[1]
		dprtype=obs[2]
		dataset_id_eso=obs[3]
		dataset_raw=obs[4]
		dataset_red=obs[5]
		print("Now processing: ", dataset_id_eso)
		XQC_plot(dataset_id_eso, dataset_raw, dataset_red, dprtype, night,arm)

def QC_all():
	nights=ascii.read("/Users/leo/Projekte/LP-BAT/XSHOOTER/observations/log/nights_unique.txt")
	for night in nights:
		QC_night(night[0])

QC_all()

dataset_id_eso="XSHOO.2015-05-13T09:56:59.151"
dataset_raw="/Volumes/astrodata2/XSHOODATA/2015-05-12/XSHOO.2015-05-13T09:56:59.151.fits"
dataset_red="/Users/leo/Projekte/LP-BAT/XSHOOTER/data/reflex_end_products/2015-11-12_pipeline_2-6-8/XSHOO.2015-05-13T09_56_59.151_tpl/LTT7987_onoff_IFU_FLUX_IFU_MERGE3D_DATA_OBJ_VIS.fits"
#XQC_plot(dataset_id_eso, dataset_raw, dataset_red, "FLUX", "2015-05-12")


dataset_id_eso="XSHOO.2015-05-13T10:51:59.198"
dataset_raw="/Volumes/astrodata2/XSHOODATA/2015-05-12/XSHOO.2015-05-13T10:51:59.198.fits"
dataset_red="/Users/leo/Projekte/LP-BAT/XSHOOTER/data/reflex_end_products/2015-11-12_pipeline_2-6-8/XSHOO.2015-05-13T10_51_59.198_tpl/Calibration_TELL_IFU_MERGE3D_DATA_OBJ_VIS.fits"
#XQC_plot(dataset_id_eso, dataset_raw, dataset_red, "TELL", "2015-05-12")
## this file is missing!

night="2015-05-12"
dprtype="TELL"
dataset_id_eso="XSHOO.2015-05-13T04:23:35.008"
dataset_raw="/Volumes/astrodata2/XSHOODATA/2015-05-12/XSHOO.2015-05-13T04:23:35.008.fits"
dataset_red="/Users/leo/Projekte/LP-BAT/XSHOOTER/data/reflex_end_products/2015-11-12_pipeline_2-6-8/XSHOO.2015-05-13T04_23_35.008_tpl/Hip076234_TELL_IFU_MERGE3D_DATA_OBJ_VIS.fits"
#XQC_plot(dataset_id_eso, dataset_raw, dataset_red, dprtype, night)