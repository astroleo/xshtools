from postprocess import flatten_ob
from molecfit_wrapper import *

## xshoo pipeline
##
## combines various procedures to produce final calibrated files for some OBs

#### missing
## merge to long UVB/VIS/NIR 1D spec // check cross-calibration in flux and w.r.t. telluric fit parameters (molecfit)
## shift to redshift = 0? / re-sample wavelength to log lambda?


#ob_name="MCG-06-30-015_1"
#ob_name="NGC6814_1"
#ob_name="NGC6814_2"
#ob_name="ESO208-G021_1"
#ob_name="NGC1079_1"
##
## continue here
#ob_name="NGC2110_1" ## VIS molecfit error: GDAS profile not found, molecfit crashes
ob_name="NGC2110_2"

##
## inactive gals
NGC0718_1
NGC1079_1
NGC1315_1
NGC1947_1
ESO208-G021_1
NGC2775_1
NGC3175_1
ESO093-G003_1
NGC3717_1
NGC3749_1
NGC4224_1
NGC4254_1
NGC4260_1
NGC5037_1
NGC5845_1
NGC5921_1
IC4653_1
NGC7727_1


#flatten_ob(ob_name)

run_molecfit(ob_name,"NIR")
molecfit_QC(ob_name,"NIR")
run_molecfit(ob_name,"VIS")
molecfit_QC(ob_name,"VIS")
#
flux_calibrate(ob_name,"NIR")
flux_calibrate(ob_name,"VIS")
flux_calibrate(ob_name,"UVB")
