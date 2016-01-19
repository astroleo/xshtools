##
## calculate atmospheric dispersion as given by Bönsch & Potulski, 1998
## (and compare to measured dispersion in X-SHOOTER data)

import numpy as np
from matplotlib import pyplot as plt


##
## symbols mostly as in the paper
##
def ref_index(lam_mu, x_CO2=0.0004, T_Celsius=20, p_Pascal=1000e2, f_Pascal=0.):
	## Eq. (6a) -- for reference conditions T_Celsius = 20, p_Pascal = 1000e2, x_CO2=0.0004
	n1 = 1 + 1e-8 * (8091.37 + 2333983/(130 - (1/lam_mu)**2) + 15518/(38.9 - (1/lam_mu)**2))
	##
	## Eq. (7) -- for differing CO_2 content
	n2 = 1 + (n1-1) * (1 + 0.5327 * (x_CO2 - 0.0004))
	##
	## Eq. (8) -- for deviations in temperature and pressure
	n3 = 1 + ((n2-1) * p_Pascal)/93214.6 * (1 + 1e-8 * (0.5953 - 0.009876 * T_Celsius) * p_Pascal)/(1 + 0.0036610 * T_Celsius)
	##
	## Eq. (9a) -- for moist air
	n4 = n3 - f_Pascal * (3.8020 - 0.0384 * (1/lam_mu)**2) * 1e-10
	
	return(n4)

##
## derive partial pressure of hydrogen, given relative humidity, temperature and pressure
## from  https://en.wikipedia.org/w/index.php?title=Relative_humidity&oldid=669773787
##
def f_H2O(p_Pascal, T_Celsius, RH):
	p_mbar = p_Pascal * 1e-2
	##
	## calculate equilibrium vapor pressure
	evp = (1.0007 + 3.46 * 1e-6 * p_mbar) * (6.1121) * np.exp(17.502 * T_Celsius / (240.97 + T_Celsius))
	##
	## calculate partial pressure of water vapor
	pp = evp * RH
	return(pp)

##
## Eq. (1) from menezes2014
##    definitions:
##    z: true zenith distance
##    zeta: observed zenith distance
##    R := z - zeta
##    n: refractive index close to the Earth's surface (assumed to be given by ref_index above)
##    
def DAR(zeta, lam_mu):
	##
	## standard values for Paranal from http://www.eso.org/gen-fac/pubs/astclim/lasilla/diffrefr.html
	T_Celsius=11.5
	p_Pascal=743e2
	RH=0.145

	f_Pascal=f_H2O(p_Pascal, T_Celsius, RH)

	return(206265 * (ref_index(lam_mu, T_Celsius=T_Celsius, p_Pascal=p_Pascal, f_Pascal=f_Pascal) - 1) * np.tan(zeta))

def dDAR(zeta_deg, lam_mu, lam_mu_ref):
	zeta = np.deg2rad(zeta_deg)
	DAR_ref = DAR(zeta,lam_mu_ref)
	return(DAR(zeta, lam_mu) - DAR_ref)

def plot_dDAR():
	lam_mu = 0.3+np.arange(2200)/1000
	lam_mu_ref = 0.55
	zeta_deg = 10*np.arange(8)
	plt.ylim([-3,5])
	##
	## mark X-SHOOTER arms
	plt.fill_betweenx(plt.ylim(),0.3,0.55,color="blue",alpha=0.3)
	plt.fill_betweenx(plt.ylim(),0.55,1.0,color="green",alpha=0.3)
	plt.fill_betweenx(plt.ylim(),1.0,2.5,color="red",alpha=0.3)

	for z in zeta_deg:
		am = 1/np.sin(np.deg2rad(90-z))
		plt.plot(lam_mu, dDAR(z, lam_mu, lam_mu_ref), label="{0} deg, airmass {1:5.2f}".format(z,am))
	
	plt.legend()
	plt.title("Differential atmospheric refraction w.r.t. {0:5.2f} micron".format(lam_mu_ref))
	plt.xlabel("Wavelength in micron")
	plt.ylabel("Differential atmospheric refraction in arcsec")

	plt.text(2.4, -2.5, "Paranal standard conditions: 11.5 deg C, 743 hPa, 14.5 % RH", fontsize=9, ha="right")
