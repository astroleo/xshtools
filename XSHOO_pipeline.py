from Xspec import Xspec
import sys

ob = sys.argv[1]

print(ob)
xs=Xspec(ob, delete_old=True)

# 2016-05-23: changed order of reduce and build_dataset -- QC plots are only available after reduction
xs.reduce_data()
xs.build_dataset_QC_page()

xs.intercalibrate()
#xs.write_spec_starlight_bc03()
#xs.write_SL_infile()
#xs.run_starlight()
