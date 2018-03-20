## prepare first iteration
import shutil
import os
import numpy as np

from matplotlib import pyplot as plt
from astropy.io import ascii

from Xspec import Xspec
from SLplotlib import get_sfh_summary
from sigma_clip_starlight import sigma_clip_starlight
from obscure_spectrum import obscure_spectrum

number_of_iterations=4

# from xshtools/plotSL_collage.py
def plot_sfh_summary(file_popvec,x,y,dx=800):
	young,intermediate,old=get_sfh_summary(file_popvec,young=1e8,old=1e9)
	
	lw=7
	plt.plot([x,x+dx*young/100.],[y,y],'b-',linewidth=lw)
	plt.plot([x+dx*young/100.,x+dx*(young+intermediate)/100.],[y,y],'g-',linewidth=lw)
	plt.plot([x+dx*(young+intermediate)/100.,x+dx],[y,y],'r-',linewidth=lw)


def plot_mask_evolution(ob_name):
	basedir=ob_name+"/"

	## default mask for this source
	default_mask=basedir+"default/mask"
	if os.path.isfile(default_mask):
		mask=ascii.read(default_mask,data_start=1)
		for j in np.arange(len(mask["col1"])):
			if mask[j]["col1"] > 9000:
				continue
			plt.plot([mask[j]["col1"],mask[j]["col2"]],[0,0],"k")
	
	## AGN mask and iterations on it
	for i in np.arange(number_of_iterations):
		iterationdir=basedir+"iteration"+str(i)+"/"
		maskfile=iterationdir+"mask"
		mask=ascii.read(maskfile,data_start=1)
		for j in np.arange(len(mask["col1"])):
			if mask[j]["col1"] > 9000:
				continue
			plt.plot([mask[j]["col1"],mask[j]["col2"]],[1+i,1+i],"k")
		##
		## assemble population vector for this iteration
		popvecfile=iterationdir+ob_name+".out.popvec"
		plot_sfh_summary(popvecfile,9000,1+i)

	plt.ylabel("Iteration")
	plt.title("Mask evolution")
	plt.savefig(basedir+"mask_evolution.pdf")
	
AV=0

def evolve_mask(ob_name,AV=AV):
	if AV != 0:
		print("------------------------------------")
		print("------------------------------------")
		print("------------------------------------")
		print("------------------------------------")
		print("WARNING: AV != 0 -- the original input file will be changed")
		print("These changes will be undone at the end of the script.")
		print("If the script crashes, however, the changed input spectra will still be there.")
		print("------------------------------------")
		print("------------------------------------")
		print("------------------------------------")
		print("------------------------------------")
		
	mask_dir="/Users/leo/STARLIGHT/masks"
	out_dir="/Users/leo/STARLIGHT/out_dir"
	spec_dir="/Users/leo/STARLIGHT/spectra/"
	out_file_base=out_dir+"/"+ob_name+".out"
	#base_dir=os.getenv("PROJECTS")+"/LLAMA/XSHOOTER/Ric_visit_2017/mask_evolution/"
	base_dir="/Users/leo/strwCloud/sharing/Ric_visit_2017/mask_evolution/"
	result_dir=base_dir+ob_name
	
	##
	## remove previous results
	if os.path.exists(result_dir):
		shutil.rmtree(result_dir)
		print("Removed old result dir: " + result_dir)

	##
	## add extra extinction to spectrum, if required
	if AV != 0:
		## move original spectrum to a safe place
		specfile=spec_dir+ob_name+".txt"
		specfile_orig=specfile+".orig"
		shutil.move(specfile,specfile_orig)
		## apply extinction to input spectrum
		obscure_spectrum(specfile_orig,AV,specfile)

	mask_name=mask_dir+"/mask"
	plot="/Users/leo/STARLIGHT/plot/"+ob_name+".pdf"
	results=[plot,out_file_base,out_file_base+".masks",out_file_base+".popvec",out_file_base+".synspec"]

	xs=Xspec(ob_name)
	#xs.write_SL_infile()

	##
	## STEP 1 -- DEFAULT MASK for this galaxy
	##
	default_mask=mask_dir+"/default/"+ob_name+".sm"
	if os.path.isfile(default_mask):
		shutil.copy(default_mask,mask_name)
		xs.run_starlight()
		default_result_dir=result_dir+"/default"
		if not os.path.exists(default_result_dir):
			os.makedirs(default_result_dir)
		for res in results:
			shutil.move(res,default_result_dir)
		shutil.move(mask_name,default_result_dir)

	##
	## STEP 2 -- DEFAULT AGN mask
	##
	default_agn_mask=mask_dir+"/Rogerio_masks/Rogerio_NGC3081_1.sm"
	shutil.copy(default_agn_mask,mask_name)
	xs.run_starlight()

	##
	## STEP 3 -- Start iteration (sigma clipping) on AGN mask
	##
	for i in np.arange(number_of_iterations):
		iteration_dir=result_dir+"/iteration"+str(i)+"/"
		if not os.path.exists(iteration_dir):
			os.makedirs(iteration_dir)
		shutil.move(mask_name,iteration_dir)
		sigma_clip_starlight(ob_name,mask_name,xs.z,iteration_dir)
		for res in results:
			shutil.move(res,iteration_dir)
		xs.run_starlight()
	os.chdir(base_dir)
	plot_mask_evolution(ob_name)
	
	if AV != 0:
		shutil.move(specfile,result_dir)
		shutil.move(specfile_orig,specfile)
		print("Moved modified input spectrum to result dir")
		print("Restore original unmodified input spectrum file")

if __name__ == "__main__":
    import sys
    evolve_mask(sys.argv[1])
