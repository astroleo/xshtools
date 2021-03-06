from matplotlib import pyplot as plt
from astropy.io import ascii
import numpy as np
import os
import glob

import pdb


##
## CHANGE LOG
##
## 2018-03-21   Added parameter SL_infile to plotSL (without it the call in XSpec:740 would crash it)
## 2017-08-01   Changed residual plot to show +/- 5 sigma and draw line at +/- 2 sigma



##
## wishlist
##    - widget button to show plots of next object
##    - read llow_norm, lupp_norm from StVc04.C11.config, instead of hardcoding it here
##

##
## synRange: restrict plot to wavelength range in which the synthesis was done

## Note: Xcfg does not exist!?

def plotSL(id, SL_infile, pdf=True, synRange=False, plotSFH=True, llow_norm=6810, lupp_norm=6870, Olsyn_ini=3800, Olsyn_fin=10000):
	os.chdir(os.getenv("HOME")+"/STARLIGHT")
	
	file_obs="spectra/"+id+".txt"
	file_syn="out_dir/"+id+".out.synspec"
	file_masks="out_dir/"+id+".out.masks"
	file_popvec="out_dir/"+id+".out.popvec"

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


	###
	## wavelength range in which the fit was done
	###
	if synRange == True:
		w_min=0
		w_max=0
		for line in open(Xcfg):
			if (w_min != 0) & (w_max != 0):
				break
			a=line.split()
			if a[1] == "[Olsyn_ini]":
				w_min=np.float(a[0])
			if a[1] == "[Olsyn_fin]":
				w_max=np.float(a[0])
	else:
		w_min=wobs[0]
		w_max=wobs[-1]

#	print("w_min and w_max are: {w_min},{w_max}".format(w_min=w_min,w_max=w_max))

	##
	## later: read also Mask file to plot masked emission lines
	#with open('Mask.ESO093.BN','r') as f:
	#	for line in f:
	#		if i == 1:
	#			nlines=line
	#		i+=1

	if plotSFH == False:
		gridsize=(4,3)
		gridpos_spec=(0,0)
		gridpos_res=(3,0)
	if plotSFH == True:
		gridsize=(4,4)
		gridpos_spec=(0,0)
		gridpos_res=(3,0)
		gridpos_sfh=(2,3)

	plt.subplot2grid(gridsize,gridpos_spec,rowspan=3,colspan=3)
	plt.minorticks_on()
	plt.tick_params(axis='both', which='major', labelsize=8)
	plt.title(id)

	lw=0.5
	plt.plot(wobs, fobs, 'k-', linewidth=lw)
	
	##
	## plot observed spectrum in yellow where flagged
#	plt.plot(wobs[np.where((flagobs != 0) & (wobs < 7200))], fobs[np.where((flagobs != 0) & (wobs < 7200))], 'y|')
#	plt.plot(wobs[np.where((flagobs != 0) & (wobs > 7200))], fobs[np.where((flagobs != 0) & (wobs > 7200))], 'y|')
#	plt.plot(wobs[np.where(flagobs != 0)], fobs[np.where(flagobs != 0)], 'y|', linewidth=lw)
	
	ylim_lo = 0.5
	ylim_hi = 3
	
	w_min=Olsyn_ini
	w_max=Olsyn_fin
	
	plt.xlim((w_min,w_max))
#	plt.ylim((ylim_lo,ylim_hi))
	plt.ylim((0,2))
	plt.ylabel(r'$F_{\lambda}$' + ' [normalized]',fontsize=10)

	##
	## overplot masked regions in red
	for m in masks:
		plt.plot(wobs[(wobs >= m['lam1']) & (wobs <= m['lam2'])], 
			fobs[(wobs >= m['lam1']) & (wobs <= m['lam2'])], 'r', linewidth=lw)


	plt.plot(wobs, nobs, 'g', linewidth=lw)
	plt.plot(wsyn, fsyn, 'b', linewidth=lw)

	##
	## residual plot
	try:
		assert(fsyn.size < fobs.size)
	except:
		print("fsyn.size is greater than fobs.size")

	## bring fsyn and fobs to same length
#	pdb.set_trace()
	fres = fobs[(wobs >= wsyn[0]) & (wobs <= wsyn[-1])] - fsyn
	## this is the noise cut to the same wavelength range as the residuals
	nobs_res = nobs[(wobs >= wsyn[0]) & (wobs <= wsyn[-1])]
	fres/=nobs_res

	plt.subplot2grid(gridsize,gridpos_res,rowspan=1,colspan=3)
	plt.minorticks_on()
	plt.tick_params(axis='both', which='major', labelsize=8)

	##
	## cut out masked regions
	ix=np.where(spec_syn['weight'] > 0)
	plt.plot(wsyn[ix], fres[ix], 'k.', markersize=1.0)
	##
	## mark points that have not been used for the fit
	##
	plt.plot(wsyn, np.full(np.shape(fres),-2,dtype=int),color='grey',linestyle="dotted")
	plt.plot(wsyn, np.full(np.shape(fres),2,dtype=int),color='grey',linestyle="dotted")

	plt.xlim((w_min,w_max))
	plt.ylim((-5,5))
	plt.xlabel(r'$\lambda$' + ' [$\AA$]',fontsize=10)
	plt.ylabel(r"Residual/$\sigma$",fontsize=10)
		
	if plotSFH == True:
		plt.subplot2grid(gridsize,gridpos_sfh,rowspan=2,colspan=1)
		plot_sfh(file_popvec)
		plt.subplots_adjust(left=0.08,bottom=0.09,right=0.96,top=0.93,wspace=0.55,hspace=0.32)


	if pdf == True:
		file_plot=os.getenv("HOME")+"/STARLIGHT/plot/"+id+".pdf"
		plt.savefig(file_plot)
		plt.tight_layout()
		plt.close()
	else:
		plt.show()

###
## get population vector as functions of light and mass
##
## nice to have:
##    - more generic (e.g. automatically determine no. of metallicities from output data)
##    - check whether age distributions are really the same for all metallicities
##    - distinguish light/mass weighted
###
def plot_sfh(file_popvec):
	popvec=np.genfromtxt(file_popvec)

	## generate popvecs (light/mass weighted) for the three metallicities
	popvec_sub=popvec[np.where(popvec[:,5]<0.02)]
	popvec_sol=popvec[np.where(popvec[:,5]==0.02)]
	popvec_sup=popvec[np.where(popvec[:,5]>0.02)]

	popvec_lw_sub=popvec_sub[:,1]
	popvec_lw_sol=popvec_sol[:,1]
	popvec_lw_sup=popvec_sup[:,1]

	agevec_log=np.log10(popvec_sol[:,4])

	if len(popvec_lw_sub)>0:
		## assume we have three metallicities
		plt.bar(agevec_log,popvec_lw_sub, width=0.15,color='b',label='sub-solar')
		plt.bar(agevec_log,popvec_lw_sol, width=0.15,color='g',bottom=popvec_lw_sub,label='solar')
		plt.bar(agevec_log,popvec_lw_sup, width=0.15,color='r',bottom=popvec_lw_sol+popvec_lw_sub,label='super-solar')
	else:
		##assume we have only solar
		plt.bar(agevec_log,popvec_lw_sol, width=0.15,color='g',label='solar')
	
	plt.xticks((6,7,8,9,10),('6','7','8','9','10'))


	plt.legend(fontsize='xx-small')

	plt.xlabel("Log age [yr]")
	plt.ylabel("x_j (light) [%]")


def get_sfh_summary(file_popvec,young=2.5e7,old=1.4e9):
	popvec=np.genfromtxt(file_popvec)
	y=np.where(popvec[:,4] <= young)
	o=np.where(popvec[:,4] > old)
	i=np.where((popvec[:,4] > young) & (popvec[:,4] <= old))
	x_y=np.sum(popvec[y,1])
	x_i=np.sum(popvec[i,1])
	x_o=np.sum(popvec[o,1])
	##
	## normalize
	x_total=0.01*(x_y+x_i+x_o)
	x_y/=x_total
	x_i/=x_total
	x_o/=x_total
	
	return x_y, x_i, x_o

def plot_popvec_hist():
	dir=os.getenv('HOME')+'/STARLIGHT/out_dir/'
	os.chdir(dir)
	
	x_y_active=[]
	x_i_active=[]
	x_o_active=[]
	x_y_control=[]
	x_i_control=[]
	x_o_control=[]

	for file_syn in glob.glob("*.synspec"):
		id=file_syn.split(".")[0]
	
		file_popvec=id+".out.popvec"
		x_y, x_i, x_o = get_sfh_summary(file_popvec)
	
		if (id == "NGC4593") | (id == "MCG514") | (id == "NGC1365") | (id == "NGC2110") | (id == "NGC2992") | (id == "NGC3081"):
			print(id + " (AGN)")
			print(round(x_y), round(x_i), round(x_o))
			x_y_active.append(x_y)
			x_i_active.append(x_i)
			x_o_active.append(x_o)
		else:
			print(id + " (control galaxy)")
			print(round(x_y), round(x_i), round(x_o))
			x_y_control.append(x_y)
			x_i_control.append(x_i)
			x_o_control.append(x_o)

	plt.subplot(311)
	plt.hist((x_y_active,x_y_control),color=('b','r'),label=('AGN','control'))
	plt.legend()
	plt.xlim(0,100)
	plt.ylim(0,7)
	plt.title('Young population')
	plt.ylabel('Number of galaxies')

	plt.subplot(312)
	plt.hist((x_i_active,x_i_control),color=('b','r'),label=('AGN','control'))
	plt.xlim(0,100)
	plt.ylim(0,7)
	plt.title('Intermediate-age population')
	plt.ylabel('Number of galaxies')

	plt.subplot(313)
	plt.hist((x_o_active,x_o_control),color=('b','r'),label=('AGN','control'))
	plt.xlim(0,100)
	plt.ylim(0,7)
	plt.title('Old population')
	plt.xlabel('Light fraction')
	plt.ylabel('Number of galaxies')

	plt.tight_layout()

	plt.savefig('../plot/SFH_compare.pdf')
	plt.close()

if __name__ == "__main__":
    import sys
    plotSL(sys.argv[1])
