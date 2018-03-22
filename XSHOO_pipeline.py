##
## summary of commands to process basic-reduced data
##

## Change log
#
#  2018-03-21   added sigma_clipping / mask generation to this list of commands
#  2016-05-23   changed order of reduce and build_dataset -- QC plots are only available after reduction

from Xspec import Xspec
import evolve_mask
import sys
import pdb

ob = sys.argv[1]
print(ob)

if len(sys.argv) == 2:
	## use 1.8" extraction
	extraction="large"
	xs=Xspec(ob, delete_old=True)
elif len(sys.argv) == 4:
	## use 0.6" extraction with given slitlets
	extraction="small"
	slitlet_UVB = int(sys.argv[2])
	slitlet_VIS = int(sys.argv[3])
	xs=Xspec(ob, delete_old=True, slitlet_UVB=slitlet_UVB, slitlet_VIS=slitlet_VIS)

xs.reduce_data()
###xs.build_dataset_QC_page()
xs.intercalibrate()
xs.write_spec_starlight_bc03()
xs.write_SL_infile()
xs.run_starlight()
evolve_mask.evolve_mask(ob)