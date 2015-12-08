## xshoo pipeline
##
## (1) telluric calibration: flux and noise spectrum
## (2) fudge factor determination / correction
## (3) extraction with DAR compensation, merge to long UVB/VIS/NIR 1D spec?
## (4) shift to redshift = 0? / re-sample wavelength to log lambda?
##

- flatten_spectrum for SCI, TEL, FLX
- run and apply molecfit to SCI + FLX
- flux calibration
- store final spectrum + combine masks