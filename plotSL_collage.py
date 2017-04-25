from matplotlib import pyplot as plt
from astropy.io import ascii
from astropy.modeling import models, fitting
import numpy as np
import os
import glob
from scipy.signal import medfilt
from SLplotlib import get_sfh_summary

import pdb

def offset_from_spec(wobs,fobs):
	l_init = models.Linear1D()
	fit_l = fitting.LinearLSQFitter()
	l = fit_l(l_init,wobs,fobs)
	s = l.slope.value
	i = l.intercept.value
	offset = -i + 1000*s
	
	return(offset)

def plot_sfh_summary(file_popvec,x,y,dx=800):
	young,intermediate,old=get_sfh_summary(file_popvec,young=1e8,old=1e9)
	
	lw=7
	plt.plot([x,x+dx*young/100.],[y,y],'b-',linewidth=lw)
	plt.plot([x+dx*young/100.,x+dx*(young+intermediate)/100.],[y,y],'g-',linewidth=lw)
	plt.plot([x+dx*(young+intermediate)/100.,x+dx],[y,y],'r-',linewidth=lw)


llow_norm=6810
lupp_norm=6870
Olsyn_ini=3800
#Olsyn_fin=10000
Olsyn_fin=12500

obs = ascii.read("/Users/leo/Desktop/inactive_plot.txt",names=("id","offset"))

for id,offset in zip(obs["id"],obs["offset"]):
	file_obs="/Users/leo/STARLIGHT/spectra/"+id+".txt"
	file_syn="/Users/leo/STARLIGHT/out_dir/"+id+".out.synspec"
	file_masks="/Users/leo/STARLIGHT/out_dir/"+id+".out.masks"
	file_popvec="/Users/leo/STARLIGHT/out_dir/"+id+".out.popvec"

	spec_obs=ascii.read(file_obs, names=('wave','flux','noise','flag'))
	spec_syn=ascii.read(file_syn, names=('wave','fobs','fsyn','weight'))
	masks=ascii.read(file_masks, names=('lam1','lam2'))
	##
	## normalize spectrum the same way STARLIGHT does
	norm=np.median(spec_obs['flux'][(spec_obs['wave'] >= llow_norm) & (spec_obs['wave'] <= lupp_norm)])
	wobs=spec_obs['wave']
	fobs=spec_obs['flux']/norm
	nobs=spec_obs['noise']/norm
	flagobs=spec_obs['flag']

	wsyn=spec_syn['wave']
	fsyn=spec_syn['fsyn']

	lw=0.5
	medfiltwidth = 11
#	offset = offset_from_spec(wobs,fobs)
	print(id,offset)
	plt.plot(wobs, offset+medfilt(fobs,medfiltwidth), 'k-', linewidth=lw, label="")

	##
	## overplot masked regions in red
	for m in masks:
		plt.plot(wobs[(wobs >= m['lam1']) & (wobs <= m['lam2'])], 
			offset+medfilt(fobs[(wobs >= m['lam1']) & (wobs <= m['lam2'])],medfiltwidth), 'r', linewidth=lw, label="")

#	labelstr = id + "(" + str(H) + ")"
	labelstr = id
	plt.plot(wsyn, offset+medfilt(fsyn,medfiltwidth), 'b-', linewidth=lw)
	id_label_y = offset+fobs[-1]
	plt.text(11050,id_label_y,labelstr,fontsize=8)
	plot_sfh_summary(file_popvec,10100,id_label_y)

w_min=Olsyn_ini
w_max=Olsyn_fin

plt.xlim((w_min,w_max))
plt.ylim((-1.5,4.5))
plt.ylabel(r'$F_{\lambda}$' + ' [normalized]')

#plt.legend(bbox_to_anchor=(1.05, 1), loc=2, fontsize=10) ## places legend outside of axes and unfortunately also outside of plotting area
#plt.legend(loc=4, fontsize=9)

plt.savefig("plot_collage_8_9.pdf")
plt.clf()