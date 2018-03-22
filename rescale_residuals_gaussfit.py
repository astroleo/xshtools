from astropy.io import ascii
from astropy.modeling import models, fitting
import numpy as np
from matplotlib import pyplot as plt
import pdb

def stddev_residuals(f_synspec,f_plot=False):
	a=ascii.read(f_synspec,names=["wave","f_obs","f_model","weight"])
	wave=a["wave"]
	f_obs=a["f_obs"]
	f_model=a["f_model"]
	weight=a["weight"]
	res_sigma=(f_obs-f_model)*weight
	
	##
	## weight=0 means we have a masked region; let's not look at residuals there
	hist=plt.hist(res_sigma[res_sigma!=0])
	n=hist[0]
	bin_edges=hist[1]
	bin_centres=np.zeros(len(n))
	for j in np.arange(len(n)):
		bin_centres[j]=(bin_edges[j]+bin_edges[j+1])/2
	
	fit_g=fitting.LevMarLSQFitter()
	bound_parameters = {'stddev': [0,2]}
	g_init=models.Gaussian1D(amplitude=1500, mean=0, stddev=0.5, bounds = bound_parameters)
	g=fit_g(g_init,bin_centres,n)
#	print(g)
		
	x=np.linspace(-3,3,200)
	plt.plot(x,g(x))
	stddev_fitted=g.stddev.value
	plt.xlim([-3,3])
	ylim=plt.ylim()
	plotstring="stddev = {0:5.2f}".format(stddev_fitted)
	plt.text(2.9,ylim[1]*0.9,plotstring,ha="right")
	plt.xlabel("Residuals/sigma_pipeline")
	plt.title(f_synspec,fontsize=9)
	
	if f_plot:
		plt.savefig(f_plot)
		plt.clf()
	
	return stddev_fitted
	