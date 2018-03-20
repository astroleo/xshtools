from astropy.io import ascii
from screen_extinction import screen_extinction as ex
##
## multiply a given spectrum with a given extinction
##
## crop_range: crop spectra to actually used spectral range
##
def obscure_spectrum(infile,AV,outfile,crop_range=[2800,24000]):
	spec=ascii.read(infile,names=(["w","f","e","m"]))
	##
	## crop to actually used spectral range
	spec_crop = spec[(spec["w"]>crop_range[0]) & (spec["w"]<crop_range[1])]
	##
	## compute wavelength-dependent extinction
	spec_crop["f"] = ex(spec_crop["w"],AV) * spec_crop["f"]
	spec_crop["e"] = ex(spec_crop["w"],AV) * spec_crop["e"]
	##
	## ...and write it all out again
	spec_crop.write(outfile, format="ascii.fixed_width_no_header", delimiter=None)
	print("Wrote modified spectrum to " + outfile)
	