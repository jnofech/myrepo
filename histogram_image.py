import matplotlib.pyplot as plt
import numpy as np
#get_ipython().magic(u'matplotlib inline')

import astropy.io.fits as fits


def hist_make():

	data = fits.getdata("m33.ico.fits")
	data = data.squeeze()

	datanoise = np.isclose(data,-2.02482640e-07)  	#This is an array of the "noise" data entries.
	data2 = data - 100*datanoise
	dataclean = data2[data2>-100]                 	#This is a 1D array of the "not-noisy" entries of "data".


	plt.hist(dataclean.ravel(),bins=200)
	plt.axis([-30,30,0,160000])			
	plt.savefig('histogram.png')

