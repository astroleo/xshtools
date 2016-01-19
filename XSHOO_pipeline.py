from postprocess import flatten_ob
from molecfit_wrapper import *
## xshoo pipeline
##

ob_name="NGC3783_2"

#flatten_ob(ob_name)

#run_molecfit(ob_name,"NIR")
#molecfit_QC(ob_name,"NIR")
#run_molecfit(ob_name,"VIS")
#molecfit_QC(ob_name,"VIS")

#flux_calibrate(ob_name,"NIR")
#flux_calibrate(ob_name,"VIS")
flux_calibrate(ob_name,"UVB")

#### missing
## merge to long UVB/VIS/NIR 1D spec // check cross-calibration in flux and w.r.t. telluric fit parameters (molecfit)
## shift to redshift = 0? / re-sample wavelength to log lambda?
