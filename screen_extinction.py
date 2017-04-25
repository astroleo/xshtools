from astropy.io import ascii
import numpy as np
from scipy.interpolate import interp1d
import os

##
## returns factor to multiply intrinsic spectrum with to get a screen extinction of a certain A_V at lam_norm
##
## INPUTS
##    lam   wavelength vector of spectrum to be reddened (to interpolate extinction for those points)
##    AV    extinction in the visual band (or at other band if lam_norm is set, too)
##
## WAVELENGTH UNITS in this function: Angstrom
##
## CHANGE LOG
##    2016-02-01   adapted from old IDL code used for B15
##
def screen_extinction(lam, AV, lam_norm=5500):
	##
	## read and interpolate dust mass extinction curve
	##
	mrn=os.getenv("XSHTOOLS") + "/mrn_kappa_ext.dat"
	a=ascii.read(mrn,data_start=3,names=(["lam_mu","kappa_ext"]))
	##
	## convert to AA and log-scale (in order to do interpolation in log-space)
	##
	log_tau_lam = np.log10(a["lam_mu"] * 1e4)
	log_tau_val = np.log10(a["kappa_ext"])
	log_lam = np.log10(lam)
	log_lam_norm = np.log10(lam_norm)
	##
	## interpolate tau for given wavelengths and at normalization wavelength
	##
	interp_fct = interp1d(log_tau_lam,log_tau_val)
	log_tau = interp_fct(log_lam)
	log_tau_lam_norm = interp_fct(log_lam_norm)
	##
	## convert back to linear scale and normalize: tau_n = 1 at lam_norm
	##
	tau_n = 10**(log_tau - log_tau_lam_norm)
	##
	## A (mag) = 2.5 * tau * alog10(exp(1)) ~ 1.09 * tau
	## <==> tau ~ 0.92 * A (mag)
	##
	e_norm = np.exp(- AV * tau_n * 1/(2.5 * np.log10(np.exp(1))))
	
	return(e_norm)

##
## test:
## w=3000+100*np.arange(100)
## lam_norm=5500
## e=screen_extinction.screen_extinction(w,1, lam_norm=lam_norm)
## e[w==lam_norm] ## should give extinction corresponding to 1 mag, i.e. 1/2.5 = 0.4