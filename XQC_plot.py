# QC plot for each obs / raw data
#- dataset name
#- hist of raw data counts
#- plot of spectrum
#- total number of bad pixels

##
## produce a Quality Check plot for an X-SHOOTER reduced data set
##
import matplotlib
import numpy as np
from astropy.io import ascii,fits
from matplotlib import pyplot as plt
import os

def XQC_plot(dataset_id_eso, dataset_raw, dataset_red, dprtype, night, arm, ix, flux_cube, qual_cube, plot_type="simple"):
	## check if required data exist
	if not os.path.isfile(dataset_raw):
		print("RAW file missing: ", dataset_raw)
		return()
	if not os.path.isfile(dataset_red):
		print("REDUCED file missing: ", dataset_red)
		return()

	## check if plot already exists
	dir_out=os.getenv("XDIR")+"/QC/plot/"+night
	if not os.path.isdir(dir_out):
		os.mkdir(dir_out)
	if plot_type == "full":
		f_out=dir_out+"/"+dataset_id_eso.replace(":","_")+"_"+dprtype+"_"+arm+".png"
	elif plot_type == "simple":
		f_out=dir_out+"/"+dataset_id_eso.replace(":","_")+"_"+dprtype+"_"+arm+"_simple.png"
	if os.path.isfile(f_out):
		print("QC PLOT already exists: ", f_out)
		return()

	## data preparation
	hdu_raw=fits.open(dataset_raw)
	d=hdu_raw[0].data
	ob_name=hdu_raw[0].header["HIERARCH ESO OBS NAME"]
	am=hdu_raw[0].header["HIERARCH ESO TEL AIRM START"]
	fwhm=hdu_raw[0].header["HIERARCH ESO TEL AMBI FWHM START "]	
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
	
	### BEGIN OF full PLOT ###
	if plot_type == "full":
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

	### BEGIN OF simple PLOT ###
	elif plot_type == "simple":
		# plot of image
		img_collapsed=np.sum(flux_cube,axis=0)
		ix_center=ix[10000,:,:]
		aperture_borders=np.where(np.diff(ix_center[:,0])!=0)
		
		plt.imshow(img_collapsed,origin="lower")
		for y in aperture_borders:
			plt.plot([0,2],[y,y], color="white", linewidth=2)
		plt.axis("off")
		
		plt.savefig(f_out)
		print("Saved ", f_out)
		plt.close()
