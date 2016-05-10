
# 5.10.16 - Working with Subcubes from "paws_norot.fits", in spatial and spectral dimensions.


from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import scipy.stats as ss
import math

print('Welcome to Cubes2D! \n \nAvailable functions: \n  cubegen - generates a subcube from each map (12CO and 13CO). \n  structgen - generates structure function maps from a given subcube. \n  mapgen - generates 2D maps of the structure functions, for 12CO and 13CO. \n  plotgen - generates 1D plots of structure function versus radius. Shading optional.')

def cubegen(ymin,ymax,xmin,xmax, deltaX=30):
	"""Generates a subcube of the specified dimensions from the .fits files,
	   for 12CO and 13CO. Returns the subcubes for 12CO and 13CO, respectively.

	   Argument format: "(ymin,ymax, xmin,xmax.)"
	   ^ These are the parameters of the desired subcubes."""

	cube12 = SpectralCube.read("paws-30m-12co10-23as-cube.fits")
	cube13 = SpectralCube.read("paws-30m-13co10-23as-cube.fits")
	subcube12 = cube12[:,ymin:ymax,xmin:xmax]
	subcube13 = cube13[:,ymin:ymax,xmin:xmax]

	return subcube12,subcube13

def structgen(subcube12,subcube13, deltaX=30):

	"""Generates a structure function map from two given subcubes; the first
	   for 12CO and the second for 13CO.
	
	   Argument format: "(<subcube12>, <subcube13>, deltaX (default=30))".
	   "deltaX" (default: 30) is the maximum value of dX and dY."""

	moment012 = subcube12.moment(0,axis=0)
	moment013 = subcube13.moment(0,axis=0)
	moment012 = subcube12.apply_numpy_function(np.nanmax,axis=0)*u.K
	moment013 = subcube13.apply_numpy_function(np.nanmax,axis=0)*u.K

	dX = deltaX                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	nmax = abs(2*dX)+1
	S_2_12 = np.zeros([nmax,nmax])
	S_2_13 = np.zeros([nmax,nmax])

	n12,m12 = moment012.shape     # Each 2D subcube has n rows and m columns.
	n13,m13 = moment013.shape


	for dx in range (-dX,dX+1):
	    for dy in range (0,dY+1):
		
		M12 = moment012.value         			# This is the array of the subcube's values (unitless).
		P12 = np.arange(n12*m12).reshape(n12,m12) 	# This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		D12 = np.zeros([n12,m12])   			# This will be the difference between M(r) and M(r+dr).
		
		RollMap12 = np.roll(np.roll(M12,-dy,axis=0),-dx,axis=1)
		D12 = M12 - RollMap12
		
		goodpix12 = (P12 - np.roll(P12,-dy,axis=0) == -dy*m12) * (P12 - np.roll(P12,-dx,axis=1) == -dx)
		        # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		        #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		        #        that of a pixel "dy" rows below is simply dy*m.
		        # In "goodpix", pixels that have wrapped around are treated as "False".
		        
		OrigMapPower12 = (np.nanmean(M12[goodpix12]**2))
		RollMapPower12 = (np.nanmean(RollMap12[goodpix12]**2))
		
		S_2_12[(dy+dY,dx+dX)] = (np.nanmean(D12[goodpix12]**2)) / (OrigMapPower12 + RollMapPower12)
		S_2_12[-(dy+dY+1),-(dx+dX+1)] = S_2_12[(dy+dY,dx+dX)]
		

		M13 = moment013.value         			# This is the array of the subcube's values (unitless).
		P13 = np.arange(n13*m13).reshape(n13,m13) 	# This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		D13 = np.zeros([n13,m13])   			# This will be the difference between M(r) and M(r+dr).
		
		RollMap13 = np.roll(np.roll(M13,-dy,axis=0),-dx,axis=1)
		D13 = M13 - RollMap13
		
		goodpix13 = (P13 - np.roll(P13,-dy,axis=0) == -dy*m13) * (P13 - np.roll(P13,-dx,axis=1) == -dx)
		        # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		        #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		        #        that of a pixel "dy" rows below is simply dy*m.
		        # In "goodpix", pixels that have wrapped around are treated as "False".
		        
		OrigMapPower13 = (np.nanmean(M13[goodpix13]**2))
		RollMapPower13 = (np.nanmean(RollMap13[goodpix13]**2))
		
		S_2_13[(dy+dY,dx+dX)] = (np.nanmean(D13[goodpix13]**2)) / (OrigMapPower13 + RollMapPower13)
		S_2_13[-(dy+dY+1),-(dx+dX+1)] = S_2_13[(dy+dY,dx+dX)]
	return S_2_12, S_2_13


def mapgen(S_2_12, S_2_13, deltaX=30):
	"""Generates and saves a 2D colour plot of the structure functions versus
	   position; one for 12CO and one for 13CO.

	   Argument format: (S_2_12, S_2_13, deltaX). Plots are created using the resulting
	   matrices from structgen, and the same deltaX that was used in structgen."""

	dX = deltaX                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	nmax = abs(2*dX)+1


	# 2D Display
	plt.figure(1)
	plt.subplot(121)
	plt.imshow(S_2_12, interpolation = 'none', extent = [-dX,dX,-dY,dY], vmin=0, vmax=max(S_2_12.max(),S_2_13.max()), aspect='auto')
	plt.title('S_2 for 12CO')
	plt.xlabel('dx')
	plt.ylabel('dy')

	plt.subplot(122)
	plt.imshow(S_2_13, interpolation = 'none', extent = [-dX,dX,-dY,dY], vmin=0, vmax=max(S_2_12.max(),S_2_13.max()), aspect='auto')
	plt.title('S_2 for 13C0')
	plt.xlabel('dx')
	plt.ylabel('dy')
	plt.colorbar()
	plt.savefig('map1213V3.png')

def plotgen(S_2_12, S_2_13, deltaX=30, SHADING=False):
	"""Generates and saves a 1D plot of the average structure function versus 
	   radius, for each of 12CO and 13CO.

	   Argument format: (S_2_12, S_2_13, deltaX, SHADING=False). Plots are 
	   created using the resulting matrices from structgen, and the same deltaX
	   that was used in structgen.

	   NOTE: If SHADING is set to "True", then a shaded plot will be generated.
	   There are some bugs, however; so False (default) is recommended."""

	dX = deltaX                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	nmax = abs(2*dX)+1

	# Goal: Create a 1D plot, for each of 12CO and 13CO, of the average value of structure function (inside a thin ring
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
	struct_funct12, edges12, counts12 = ss.binned_statistic(
	    radiusmap[radiusmap<maxradius], S_2_12[radiusmap<maxradius], statistic=np.nanmean, bins = reselements)
	std12, edges12, counts12 = ss.binned_statistic(
	    radiusmap[radiusmap<maxradius], S_2_12[radiusmap<maxradius], statistic=np.std, bins = reselements)

	struct_funct13, edges13, counts13 = ss.binned_statistic(
	    radiusmap[radiusmap<maxradius], S_2_13[radiusmap<maxradius], statistic=np.nanmean, bins = reselements)
	std13, edges13, counts13 = ss.binned_statistic(
	    radiusmap[radiusmap<maxradius], S_2_13[radiusmap<maxradius], statistic=np.std, bins = reselements)

		# PLOTTING
	# No shading
	if SHADING==False:
		plt.figure(2)
		plt.plot(np.arange(reselements)/mult,struct_funct12,'r.',label='12CO')
		plt.plot(np.arange(reselements)/mult,struct_funct12+std12,'r:')
		plt.plot(np.arange(reselements)/mult,struct_funct12-std12,'r:')
		plt.plot(np.arange(reselements)/mult,struct_funct13,'b--',label='13CO')
		plt.plot(np.arange(reselements)/mult,struct_funct13+std13,'b-')
		plt.plot(np.arange(reselements)/mult,struct_funct13-std13,'b-')
		plt.title('Average Structure Function vs. Radial "Distance" from Center of S_2 Plots')
		plt.xlabel(' "Radius" ')
		plt.ylabel('Average S_2')
		plt.legend(loc='upper left')
		plt.savefig('plot1213_noshadingV3.png')
	else:
	# Yes, shading
		plt.figure(3)
		plt.plot(np.arange(reselements)/mult,struct_funct12,'k.',label='12CO average')
		plt.plot(np.arange(reselements)/mult,struct_funct13,'k:',label='13CO average')

		plt.fill_between(np.arange(reselements)/mult, struct_funct12-std12, struct_funct13+std13, where=struct_funct13+std13 >= struct_funct12-std12, facecolor='purple')
		plt.fill_between(np.arange(reselements)/mult, struct_funct12-std12, struct_funct13+std13, where=struct_funct12+std12 >= struct_funct13-std13, facecolor='purple')

		plt.fill_between(np.arange(reselements)/mult, struct_funct12+std12, struct_funct13+std13, where=struct_funct12+std12 >= struct_funct13+std13, facecolor='red')
		plt.fill_between(np.arange(reselements)/mult, struct_funct12-std12, struct_funct13-std13, where=struct_funct12-std12 <= struct_funct13-std13, facecolor='red')

		plt.fill_between(np.arange(reselements)/mult, struct_funct12+std12, struct_funct13+std13, where=struct_funct13+std13 >= struct_funct12+std12, facecolor='blue')
		plt.fill_between(np.arange(reselements)/mult, struct_funct12-std12, struct_funct13-std13, where=struct_funct13-std13 <= struct_funct12-std12, facecolor='blue')

		plt.fill_between(np.arange(reselements)/mult, struct_funct12-std12, struct_funct13+std13, where=struct_funct13+std13 <= struct_funct12-std12, facecolor='white')
		plt.fill_between(np.arange(reselements)/mult, struct_funct13-std13, struct_funct12+std13, where=struct_funct13-std13 >= struct_funct12+std12, facecolor='white')

		# The region within the 12CO average's standard deviation is RED. The region within the 13CO average's standard deviation is BLUE.
		# The overlapping region is PURPLE.

		plt.title('Average Structure Function vs. Radial "Distance" from Center of S_2 Plots')
		plt.xlabel(' "Radius" ')
		plt.ylabel('Average S_2')
		plt.legend(loc='upper left')
		plt.savefig('plot1213_yesshadingV3.png')
