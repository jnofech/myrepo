
# 5.10.16 - Working with Subcubes from "paws_norot.fits", in spatial dimensions only.
 

from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import scipy.stats as ss
import math

print('Welcome to Cubes2D! \n \nAvailable functions: \n  cubegen - generates a subcube from the selected map. \n  structgen - generates structure function maps from a given subcube. \n  mapgen - generates a 2D map of the structure function. \n  plotgen - generates 1D plots of structure function versus radius. Shading optional.')

def cubegen(ymin,ymax,xmin,xmax, deltaX=30):
	"""Generates a subcube of the specified dimensions from the specified
	   .fits file.

	   Argument format: "(ymin,ymax, xmin,xmax.)"
	   ^ These are the parameters of the desired subcube."""

	cube = SpectralCube.read("paws-30m-12co10-23as-cube.fits")
	subcube = cube[:,ymin:ymax,xmin:xmax]

	return subcube

def structgen(subcube, deltaX=30):

	"""Generates a structure function map from a given subcube.
	
	   Argument format: "(<subcube>, deltaX (default=30))".
	   "deltaX" (default: 30) is the maximum value of dX and dY."""

	moment0 = subcube.moment(0,axis=0)
	moment0 = subcube.apply_numpy_function(np.nanmax,axis=0)*u.K

	dX = deltaX                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	nmax = abs(2*dX)+1
	S_2 = np.zeros([nmax,nmax])

	n,m = moment0.shape     # Each 2D subcube has n rows and m columns.


	for dx in range (-dX,dX+1):
	    for dy in range (0,dY+1):
		
		M = moment0.value         			# This is the array of the subcube's values (unitless).
		P = np.arange(n*m).reshape(n,m) 	# This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		D = np.zeros([n,m])   			# This will be the difference between M(r) and M(r+dr).
		
		RollMap = np.roll(np.roll(M,-dy,axis=0),-dx,axis=1)
		D = M - RollMap
		
		goodpix = (P - np.roll(P,-dy,axis=0) == -dy*m) * (P - np.roll(P,-dx,axis=1) == -dx)
		        # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		        #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		        #        that of a pixel "dy" rows below is simply dy*m.
		        # In "goodpix", pixels that have wrapped around are treated as "False".
		        
		OrigMapPower = (np.nanmean(M[goodpix]**2))
		RollMapPower = (np.nanmean(RollMap[goodpix]**2))
		
		S_2[(dy+dY,dx+dX)] = (np.nanmean(D[goodpix]**2)) / (OrigMapPower + RollMapPower)
		S_2[-(dy+dY+1),-(dx+dX+1)] = S_2[(dy+dY,dx+dX)]
		
	return S_2


def mapgen(S_2, deltaX=30):
	"""Generates and saves a 2D colour plot of the structure functions versus
	   position.

	   Argument format: (S_2, deltaX). Plot is created using the resulting
	   matrices from structgen, requiring the same deltaX used in structgen."""

	dX = deltaX                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	nmax = abs(2*dX)+1


	# 2D Display
	plt.figure(1)
	plt.imshow(S_2, interpolation = 'none', extent = [-dX,dX,-dY,dY], vmin=0, vmax=S_2.max(), aspect='auto')
	plt.title('S_2 for selected map')
	plt.xlabel('dx')
	plt.ylabel('dy')
	plt.savefig('map2D.png')

def plotgen(S_2, deltaX=30, SHADING=False):
	"""Generates and saves a 1D plot of the average structure function versus 
	   radius.

	   Argument format: (S_2, deltaX, SHADING=False). Plots are 
	   created using the resulting matrices from structgen, requiring the same
	   deltaX that was used in structgen.

	   NOTE: If SHADING is set to "True", then a shaded plot will be generated.
	   There may be some bugs, however; so False (default) is recommended."""

	dX = deltaX                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	nmax = abs(2*dX)+1

	# Goal: Create a 1D plot of the average value of structure function (inside a thin ring
	#       at radius r) versus radius. One plot includes a shaded region indicating standard deviation.
	x = np.linspace(-dX,dX,nmax)
	y = np.linspace(-dY,dY,nmax)
	xx, yy = np.meshgrid(x,y)

	maxradius = math.floor( (dX**2 + dY**2)**0.5 )
	mult = 1                        # Increases or decreases the numbers of bins used. Most accurate results at mult=1.
	reselements = math.floor(mult*maxradius)
		                        # This is the number of "resolution elements" (i.e. the number of points
		                        #      on the struct_funct vs. radius plot) that we're dealing with.

	radiusmap = (xx**2 + yy**2)**0.5
	struct_funct, edges, counts = ss.binned_statistic(
	    radiusmap[radiusmap<maxradius], S_2[radiusmap<maxradius], statistic=np.nanmean, bins = reselements)
	std, edges, counts = ss.binned_statistic(
	    radiusmap[radiusmap<maxradius], S_2[radiusmap<maxradius], statistic=np.std, bins = reselements)

		# PLOTTING
	# No shading
	if SHADING==False:
		plt.figure(2)
		plt.plot(np.arange(reselements)/mult,struct_funct,'r.',label='CO')
		plt.plot(np.arange(reselements)/mult,struct_funct+std,'r:')
		plt.plot(np.arange(reselements)/mult,struct_funct-std,'r:')
		plt.title('Average Structure Function vs. Radial "Distance" from Center of S_2 Plots')
		plt.xlabel(' "Radius" ')
		plt.ylabel('Average S_2')
		plt.legend(loc='upper left')
		plt.savefig('plot2D_noshading.png')
	else:
	# Yes, shading
		plt.figure(3)
		plt.plot(np.arange(reselements)/mult,struct_funct,'k.',label='S_2 Average')

		plt.fill_between(np.arange(reselements)/mult, struct_funct+std, struct_funct-std, facecolor='pink')

		plt.title('Average Structure Function vs. Radial "Distance" from Center of S_2 Plots')
		plt.xlabel(' "Radius" ')
		plt.ylabel('Average S_2')
		plt.legend(loc='upper left')
		plt.savefig('plot2D_yesshading.png')
