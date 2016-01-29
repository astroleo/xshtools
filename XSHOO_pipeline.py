from postprocess import flatten_ob
from molecfit_wrapper import *
from datasets import dataset

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
##ob_list=[
###	"NGC0718_1", -- problems
#	"NGC1079_1",
#	"NGC1315_1",
###	"NGC1947_2", -- problems
#	"NGC2775_1", -- not reduced?
#	"ESO208-G021_1",
#	"NGC3175_1",
#	"NGC3351_1", # -- empty spectrum??
###	"ESO093-G003_1", ## -- molecfit problem
#	"NGC3717_2",
#	"NGC3749_1",
#	"NGC4224_1",
#	"NGC5037_1",
###	"NGC5845_1", -- problems
#	"NGC5921_1",
#	"IC4653_1" -- problems
#	"NGC7727_1" -- problems?
##	]


##
## Seyfert 2
ob_list=[
#	"NGC2110_1", -- done
#	"NGC2110_2", -- molecfit stuff missing?
#	"NGC3081_1", -- done
#	"NGC3081_2", -- done
#	NGC4388 -- not yet observed
##	"NGC5128_1" ## -- double check that extraction happens at NIR position!! / no FLX taken!! / NGC5128_2 not observed yet
#	ESO021-G004 -- not yet observed
#	NGC5506 -- same night as NGC5128_1 / no FLX observed, do reduction manually
#	NGC5728_1 ## no NIR observation?
#	"NGC5728_2", -- problem with NIR part of molecfit
#	"ESO137_1", -- done
#	"ESO137_2", -- not good according to ESO?
#	"ESO137_3", -- problems with molecfit/NIR
#	"NGC7172_1", ## -- data OK? (has been re-observed but in same night as 7172_2 -- no FLUX observed) /// reduction fails / NIR molecfit
#	NGC7172_2 -- no FLUX observed
#	NGC7172_3 -- only telluric of the night is heavily saturated
#	NGC7582 -- not reduced
]

# problematic OBs: 	"NGC1947_1": no FLUX -> use _2 instead
# 	"NGC2775_1", --> no data / observed 2015-11 / check!!
# ngc3351_1 (and _2 -- both taken in the same night): flux standard has some high counts
# 3717_1: very badly centered
# 4224_1: nearest telluric saturated, take telluric 3 hours after obs
# 4254 / 4260: not yet observed
# 7727: only telluric over-exposed

for ob_name in ob_list:
	print(ob_name)
	d=dataset(ob_name)
	d.make_QC_page()

	flatten_ob(ob_name)

	run_molecfit(ob_name,"NIR")
	molecfit_QC(ob_name,"NIR")
	run_molecfit(ob_name,"VIS")
	molecfit_QC(ob_name,"VIS")
	##
	flux_calibrate(ob_name,"NIR")
	flux_calibrate(ob_name,"VIS")
	flux_calibrate(ob_name,"UVB")
