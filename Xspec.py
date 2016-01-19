import os

class Xspec(object)
	"""
	A class to handle X-SHOOTER spectra
	
	Among other things, this does
	- build file list for a given OB (depends on config files)
	- extracts spectra from pipeline reduced data with a given aperture and taking into account differential atmospheric refraction
	- flux calibration
	- telluric calibration
	- merging of the different arms
	- for each step, "quality control" plots can be produced in order to check what happens
	"""
	def __init__(self, OB_name):
		xdir=os.getenv("PROJECTS")+"/LP-BAT/XSHOOTER/"
		base_dir_red = xdir+'data/reflex_end_products/2015-11-12_pipeline_2-6-8/'
		arms=["UVB","VIS","NIR"]

		## check if data is there
		for arm in arms:
			dataset_definition = xdir+'data/dataset_definition/'+ob_name+'_'+arm+'.txt'
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
				dir_cube=base_dir_red+dpid.replace(":","_")+"_tpl/"
				f_cube=dir_cube+ob+"_"+dpr+"_IFU_MERGE3D_DATA_OBJ_"+arm+".fits"

				if not os.path.isfile(f_cube):
					print("Cube file does not exist at", f_cube)
					continue

				dir_out=xdir+"data/spectra/"+ob_name+"/"
				if not os.path.isdir(dir_out):
					os.mkdir(dir_out)
		
				f_out=dir_out+dpr+"_"+arm+"_spec.fits"
				if os.path.isfile(f_out):
					print("Outfile (spectrum) exists:",f_out)
					continue
	
	
	
	def get_file(dpid):
	dir=os.getenv("PROJECTS")+"/LP-BAT/XSHOOTER/data/reflex_end_products/2015-11-12_pipeline_2-6-8/"
	subdir=dpid.replace(":","_")+"_tpl"
	f=glob.glob(dir+subdir+"/*_IFU_MERGE3D_DATA_OBJ_*.fits")
	if len(f) != 1:
		raise ValueError("There is no (or more than one) merged file in sub-directory " + subdir)
	
	return(f[0])