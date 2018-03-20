from astropy.io import ascii
import numpy as np
import subprocess
import pdb
import os, shutil

from rescale_residuals_gaussfit import stddev_residuals

## Python script to perform sigma clipping for STARLIGHT results, i.e. convert STARLIGHT fit residuals into new mask

## take fit result fitted with initial mask (Rogerio_NGC3081_1.sm)
## go through residuals

# residuals: normalised (to local sigma) list of residuals
# lambdas: list of lambdas

##
## CHANGE LOG
##
## 2017-08-02  Rescale input error distribution so that ~68% of residuals are within 1 stddev.
## 2017-08-01  exclude CaT from the masking process -- it shall always be fitted
## 2017-08-01  add atmospheric A band as region that is always masked
## 2017-08-01  add option to increase size of mask by some pixels on either side

##
## read relevant data
def sigma_clip_starlight(ob_name,mask_name,redshift_source,iteration_dir):
	f_synspec="/Users/leo/STARLIGHT/out_dir/"+ob_name+".out.synspec"
	f_inspec="/Users/leo/STARLIGHT/spectra/"+ob_name+".txt"	
	synspec=ascii.read(f_synspec)
	inspec=ascii.read(f_inspec)
	## ***** ERROR RESCALING HAPPENS HERE *****
	## first find re-scaling factor for input uncertainties
	##    i.e. find out if weights (= 1/sigma_pipeline) match with 
	##    distribution of residuals
	f_plot=iteration_dir+ob_name+"_residual_distribution.png"
	stddev_res = stddev_residuals(f_synspec,f_plot=f_plot)	
	##
	## Uncertainty rescaling part 2:
	##    do the actual re-scaling of the errors, i.e. modify inspec
	inspec["col3"]*= stddev_res
	##
	## now let's make a backup of the original input spectrum
	##    if the file doesn't exist, otherwise we would be overwriting
	##    the backup with a modified file if we run this iteratively
	backupfile=f_inspec+".orig"
	if not os.path.isfile(backupfile):
		shutil.copy(f_inspec,backupfile)
	##
	## now let's overwrite the original file with the modified error spectrum
	inspec.write(f_inspec,format="ascii.no_header",overwrite="True")

	##
	## get normalisation factor fobs_norm defined as:
	## fobs/fobs_norm = f_obs_starlight
	outfile="/Users/leo/STARLIGHT/out_dir/"+ob_name+".out"
	a=subprocess.run("grep fobs_norm " + outfile + " | awk '{print $1}'",shell=True,universal_newlines=True,stdout=subprocess.PIPE)
	fobs_norm=np.float(a.stdout.strip())

	lambda_obs=inspec["col1"]
	uncertainty_normalised=inspec["col3"]/fobs_norm

	lambdas=synspec["col1"]
	flux_obs=synspec["col2"]
	flux_model=synspec["col3"]

	##
	## crop uncertainty vector to length of output vectors
	uncertainty_normalised_cropped=uncertainty_normalised[((lambda_obs>=lambdas[0]) & (lambda_obs<=lambdas[-1]))]

	residual_sigma=(flux_obs-flux_model)/uncertainty_normalised_cropped

	masks=[]
	masked=False
		
	A_band=np.float64([7570,7710]) ## define atmospheric A band region which we always want to mask
	A_band/=(1.+redshift_source)  ## A band as shifted to the restframe wavelength of the galaxy
	
	CaT=[8498,8542,8662] ## central wavelengths from Andretta et al. (2005), A&A, 430, 669
	CaT_buffer=15 ## intrinsically the half-width (90% of the feature) of the CaT lines is ca. 3 Angstrom, make it a bit wider to capture broadened feature, but not too wide as not to cover the [Fe II] line in between the 2nd/3rd line

	j=0
	##
	## need while loop here in order to be able to jump over values (e.g. when forcing a masked or non-masked region)
	while True:
		if j>=np.size(residual_sigma):
			break
		## we always want to mask the atmospheric A band (ozone band)
		## we want to catch this condition exactly once when the loop hits the beginning of the A band (which it will not do exactly because the boundaries of this band are np.float64 due to the redshift of the source, therefore the round clause)
		if lambdas[j]==np.round(A_band[0]):
			print("Now we are entering the A band")
			if not masked:
				mask_start=A_band[0]
			mask_end=A_band[1]
			masks.append([mask_start,mask_end])
			print("we have now added this mask for a source of redshift {0}".format(redshift_source))
			print(masks[-1])
			j+=np.int(np.round(np.diff(A_band)[0]))
			masked=False
			continue
			
		##
		## rational for sigma clipping: we want to find and mask emission lines which should be
		##    highly significant outliers. We then extend the mask until the residual flux is less than 1 sigma.
		##    This way we capture most of the wings of the line without making the masks too wide
		##
		## according to test (LB, 1 Aug 2017) it does not matter to STARLIGHT whether there are overlapping
		##    intervals defined in the mask file. Even multiple masks / overlaps do not matter.
		##    Masks that go beyond the data range are simply ignored.
		##
		if np.abs(residual_sigma[j]) > 5:
			## assume we have found an un-masked emission line; now follow the 
			##    residuals before that wavelength until they are < 1 sigma
			j_outlier=j
			while residual_sigma[j] > 1:
				if j<0:
					break
				j-=1
			mask_start=lambdas[j]

			## now do the same for the red wing of the line
			j=j_outlier
			#pdb.set_trace()
			while residual_sigma[j] > 1:
				if j>np.size(residual_sigma)-2:
					break
				j+=1
			mask_end=lambdas[j]
			
			masks.append([mask_start,mask_end])
		j+=1

	with open(mask_name,"w") as f:
		f.write(str(len(masks))+"\n")
		for i in masks:
			txt=str("{0} {1} 0\n".format(i[0],i[1]))
			f.write(txt)
