import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy.ndimage import gaussian_filter1d
from astropy.io import fits
import os

##
## wishlist
##    - set residual flux to 0 where spectrum is masked or flagged
##    - elegant table constructor that reads the various txt files + joins them 
##         to a big relational table to be easily accessed by all scripts; to 
##         replace read_config
##


##
## units
##    wave   Angstrom
##    flux   erg/(s cm^2 Angstrom)
##

## following APOGEE technical note by C. Allende Prieto April 8, 2011
## using parameters from Peck & Reeder 1972 (Tab. 1)
##
## INPUT
##    lam_vac   wavelength in vacuum in Angstrom
## 
## RETURNS
##    lam_air   wavelength in air in Angstrom
##
def vactoair(lam_vac):
	lam_vac_mu=lam_vac/10000.
	a=0
	b1=5.791817e-2
	b2=1.67909e-3
	c1=238.0185
	c2=57.362
	x = a + b1/(c1 - 1/lam_vac_mu**2) + b2/(c2 - 1/lam_vac_mu**2)
	lam_air_mu = lam_vac_mu / (x+1)
	lam_air = lam_air_mu * 10000.
	return lam_air

##
## get ASCII spectrum from David
##
def get_ascii_spec(f,z):
	spec = np.genfromtxt(f, dtype={'names': ('wave', 'flux', 'noise','mask'),
		'formats': ('f','f','f','f')},comments='#')
	
	## check that mask is either 0 or 2
	ix = np.where((spec['mask'] != 0) & (spec['mask'] != 2))
	try:
		assert(np.size(ix) == 0)
	except:
		print("masks are not only 0 or 2!")

	## shift to restframe wavelength
	spec['wave'] = spec['wave']/(1+z)

	return spec

##
## mask a non-redshifted (telluric) region of a de-redshifted spectrum
##
def add_mask_telluric(wave, imask, region, z):
	blue_edge = region[0]/(1+z)
	red_edge = region[1]/(1+z)
	imask[(wave > blue_edge) & (wave < red_edge)] = 2

def write_spec_starlight(wave,flux,noise,mask,fname,z):
	f_factor=1e16
	lamgrid=4000+5*np.arange(3000)

	##
	## convolve to resolution of spectral base
	##
	flux=gaussian_filter1d(flux,7)
	## find out right smoothening kernel

	##
	## interpolation is flux-conserving when accounting for different spectral pixel size
	##    dl is 0.2 in UVB and VIS, 0.6 in NIR; here we re-sample to 2 \AA
	##
	f_interpol=5.
	f=interpolate.interp1d(wave,flux)
	iflux=f_interpol*f_factor*f(lamgrid)
	
	##
	## interpolate inverse variance (since this is the additive quantity)
	noise*=f_factor
	v=interpolate.interp1d(wave,1/(noise**2))
	ivariance=v(lamgrid)
	inoise=1/np.sqrt(ivariance)
	###
	## STARLIGHT seems to have problems with very tiny errors we have at the moment (2015-05-15)
	## WORKAROUND: multiply noise by some factor until noise levels become reasonable
	###
#	inoise*=30
	
	##
	## Careful when interpolating mask
	m=interpolate.interp1d(wave,mask)
	imask=m(lamgrid)
	imask[imask > 0] = 2
	
	##
	## masks for atmospheric A and B bands (and de-redshift them)
	## from https://www.cfa.harvard.edu/~cnowlan/papers/nowlan_etal_maestro_pt_model_jqsrt2007.pdf
	##    strong A band centered at 7620 \AA
	##    weaker B band centered at 6900 \AA
	add_mask_telluric(lamgrid, imask, [7500,7700], z)
	add_mask_telluric(lamgrid, imask, [6800,7000], z)
	##
	## masks for transition regions of XSHOOTER spectral arms
	add_mask_telluric(lamgrid, imask, [5530,5600], z)
	add_mask_telluric(lamgrid, imask, [10100,10200], z)
	
	##
	## add masks for sky lines in near-IR
	## use list from Ric and mask all regions where OH line flux > 0.05
	file_oh='/Users/leo/SINFO-Tools/skysub/ohspecHK.fits'
	hdu=fits.open(file_oh)
	hdr=hdu[0].header
	waveOH_mu=(np.arange(hdr['NAXIS1'])+1-hdr['CRPIX1'])*hdr['CDELT1']+hdr['CRVAL1']
	waveOH=waveOH_mu*10000
	ohspec=hdu[0].data
	##
	## truncate to J/H band
	## TODO!! put this in mask file
	waveOH=waveOH[waveOH < 18000]
	ohspec=ohspec[waveOH < 18000]
	##
	## interpolate OH spec to lamgrid resolution and mask regions affected by OH lines
	fOH=interpolate.interp1d(waveOH,ohspec,bounds_error=False,fill_value=0)
	iohspec=fOH(lamgrid)
	imask[np.where(iohspec > 0.05)] = 2
	
	##
	## mask out entire NIR
	#imask[np.where(lamgrid > 9000)] = 2

	outdir=os.getenv('PROJECTS') + '/LP-BAT/XSHOOTER/STARLIGHT/STARLIGHTv04/spectra/'
	np.savetxt(outdir+fname,np.transpose((lamgrid,iflux,inoise,imask)),fmt=('%5.0f.   %8.3f   %8.3f    %3.0f'))
	print("Wrote {fname}".format(fname=fname))	

##
## read one or all XSHOOTER config information
##
def Xconfig(id=None):
	f_config='/Users/leo/Projekte/LP-BAT/XSHOOTER/spec_plot/specplot.cfg'
	cfgs = np.genfromtxt(f_config, 
		dtype={'names': ('id', 'OB', 'z', 'aori', 'factor_uvb', 'factor_nir'),
		'formats': ('S7','f1','f5','S1','f3','f3')},comments='#')
	if id != None:
		ix=np.where(cfgs['id'] == b''+id.encode())
		try:
			assert(np.size(ix) == 1)
		except:
			print("Object id {id} not known or not unique.".format(id=id))
		ix=ix[0][0]
		cfgs=cfgs[ix]
	return cfgs

def Xspec(id,OB):
	cfg=Xconfig(id=id)
	factor_uvb = cfg['factor_uvb']
	factor_nir = cfg['factor_nir']
	z = cfg['z']
	
	dir='/Users/leo/Projekte/LP-BAT/XSHOOTER/ascii_spectra/'
	OB=str(OB)
	f_uvb=dir+id+'_uvb_'+OB+'_integspec.ascii'
	f_vis=dir+id+'_vis_'+OB+'_integspec.ascii'
	f_nir=dir+id+'_nir_'+OB+'_integspec.ascii'

	spec_uvb = get_ascii_spec(f_uvb,z)
	wave_uvb=spec_uvb['wave']
	flux_uvb=spec_uvb['flux']
	noise_uvb=spec_uvb['noise']
	mask_uvb=spec_uvb['mask']
	w_uvb=wave_uvb[wave_uvb > 3200]
	f_uvb=flux_uvb[wave_uvb > 3200]
	n_uvb=noise_uvb[wave_uvb > 3200]
	m_uvb=mask_uvb[wave_uvb > 3200]

	spec_vis = get_ascii_spec(f_vis,z)
	wave_vis=spec_vis['wave']
	flux_vis=spec_vis['flux']
	noise_vis=spec_vis['noise']
	mask_vis=spec_vis['mask']
	w_vis=wave_vis[(wave_vis > 5530) & (wave_vis < 10200)]
	f_vis=flux_vis[(wave_vis > 5530) & (wave_vis < 10200)]
	n_vis=noise_vis[(wave_vis > 5530) & (wave_vis < 10200)]
	m_vis=mask_vis[(wave_vis > 5530) & (wave_vis < 10200)]

	spec_nir = get_ascii_spec(f_nir,z)
	wave_nir=spec_nir['wave']
	flux_nir=spec_nir['flux']
	noise_nir=spec_nir['noise']
	mask_nir=spec_nir['mask']
	w_nir = wave_nir[wave_nir > 9800]
	f_nir = flux_nir[wave_nir > 9800]
	n_nir = noise_nir[wave_nir > 9800]
	m_nir = mask_nir[wave_nir > 9800]

	wave = np.hstack((w_uvb,w_vis,w_nir))
	flux = np.hstack((factor_uvb*f_uvb,f_vis,factor_nir*f_nir))
	noise = np.hstack((factor_uvb*n_uvb,n_vis,factor_nir*n_nir))
	mask = np.hstack((factor_uvb*m_uvb,m_vis,factor_nir*m_nir))

	return wave, flux, noise, mask


def stitch_spectra(factor_uvb,factor_nir,id,z,OB):
	print("Legacy function -- will disappear soon!")
	return Xspec(id,OB)


def add_lines(l_lo, l_hi, wave, flux, details=False):
	f_linelist = '/Users/leo/Projekte/LP-BAT/XSHOOTER/spec_plot/linelist.txt'
	lines = np.genfromtxt(f_linelist,
		dtype={'names': ('id', 'l0', 'aorv', 'ref','kind'),
		'formats': ('S10','f4','S1','f4','S1')},
		comments='#')
	##
	## do vaccuum to air corrections
	i=0
	for line in lines:
		if line['aorv'].decode() == "v":
			lines['l0'][i]=vactoair(line['l0'])
			lines['aorv'][i]=b'a'
		i+=1
	##
	## sort lines by wavelength
	lines.sort(order='l0')

	xr=plt.xlim()[1]-plt.xlim()[0]
	yr=plt.ylim()[1]-plt.ylim()[0]
	ylabel_default_e = plt.ylim()[1] - 0.1 * yr
	ylabel_e = ylabel_default_e
	ylabel_default_a = plt.ylim()[0] + 0.1 * yr
	ylabel_a = ylabel_default_a
	wp=0
	
	for line in lines[(lines['l0'] > l_lo) & (lines['l0'] < l_hi)]:
		w=line['l0']
		kind=line['kind'].decode()
#		print("Line of kind {kind} at {w} Angstrom: {id}".format(kind=kind,w=w,id=id))

		dw_avg = 1
		f = np.average(flux[(wave > w - dw_avg) & (wave < w + dw_avg)])
		
		dw = w - wp
		if dw < 0.05 * xr:
			ylabel_e -= 0.1 * yr
			ylabel_a += 0.1 * yr
			if ylabel_e < 1.15 * f:
				ylabel_e = ylabel_default_e
			if ylabel_a > 0.85 * f:
				ylabel_a = ylabel_default_a
		else:
			ylabel_e = ylabel_default_e
			ylabel_a = ylabel_default_a
		
		if kind == "a":
			plt.annotate(
				line['id'].decode(),
				xy=(w, 0.9*f),
				xytext=(w, ylabel_a),
				arrowprops=dict(arrowstyle='-', linewidth=0.3, color='b', relpos=(0,0)),
				ha='left', size=5, color='b')
		else:
			plt.annotate(
				line['id'].decode(),
				xy=(w, 1.1*f),
				xytext=(w, ylabel_e),
				arrowprops=dict(arrowstyle='-', linewidth=0.3, color='b', relpos=(0,0)),
				ha='left', size=3, color='b')
		
		wp = w
