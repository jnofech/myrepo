
# 5.09.16 - Structure Function with Full Cube Information; with Command-Line Arguments


from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import scipy.stats as ss
import math


def image_make(vmin,vmax,ymin,ymax,xmin,xmax, deltaX=30, deltaV=3):
	"""Argument format: "(vmin,vmax, ymin,ymax, xmin,xmax, deltaX, deltaV.)"
		^ The first six arguments are the parameters of the desired subcube.
		WARNING: Selecting too large of a subcube will hugely increase processing time.

		"deltaX" (default: 30) is the maximum value of dX and dY. "deltaV" (default: 3)
		is the maximum value of dV.
		These latter two arguments are optional."""
	cube = SpectralCube.read('paws_norot.fits')
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	subcube = cube[vmin:vmax,ymin:ymax,xmin:xmax]

	# 1. Extracting a 3D subcube
	# 2. Compute a `data0` map (add up in the spectral direction using `data0 = subcube.filled_data[:]` Remember to take `data0.value` for what follows.
	# 3. Calculate the structure function for a small number of offsets $\delta x =\{-1,0,1\}$ and $\delta y = \{-1,0,1\}$ and $\delta z = \{-1,0,1\}$.  Given a map $M(x,y,z)$
	#
	# $$ S_2(\delta x, \delta y) = \mathrm{mean}([M(r) - M(r+\delta r)]^2)$$
	# WARNING: Selecting too large of a subcube will massively increase processing time.

	data0 = subcube.filled_data[:]

	dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dZ = deltaV		     	# Same as above, but for "dv". We call it "dZ" here since that may make it easier to visualize.
	nmax = abs(2*dX)+1
	S_2 = np.zeros([abs(2*dZ)+1,nmax,nmax])

	p,n,m = data0.shape     # The matrix "M" referred to above has n rows, m columns, and p "height units".

	# OPTIMIZED VERSION: Do NOT use. This introduces perfect rotational symmetry (wrt dx and dy) for dV!=0, but the nature of
	#    the program prevents this perfect symmetry when dV is not zero. Thus, for each dV "sheet", every value must be calculated
	#    individually; no shortcuts.

	for dz in range (0,dZ+1):
	    for dx in range (-dX,dX+1):
		for dy in range (-dY,dY+1): # OPTIMIZED version: go from 0 to dY+1.
		
		    M = data0.value         # This is the matrix "M" referred to above.
		    P = np.arange(p*n*m).reshape(p,n,m) # This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		    D = np.zeros([p,n,m])   # This will be the difference between M(r) and M(r+dr).

		    RollMap = np.roll((np.roll(np.roll(M,-dz,axis=0),-dy,axis=1)),-dx,axis=2)
		    D = M - RollMap

		    goodpix = (P - np.roll(P,-dz,axis=0) == -dz*n*m) * (P - np.roll(P,-dy,axis=1) == -dy*m) * (P - np.roll(P,-dx,axis=2) == -dx)
		            # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		            #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		            #        that of a pixel "dy" rows below is simply dy*m.
		            # In "goodpix", pixels that have wrapped around are treated as "False".
		    #goodpix = goodpix.astype('float')
		    #goodpix[goodpix==0] = 'nan'     # Now, we can disregard the wraparound pixels entirely.

		    OrigMapPower = np.nanmean(M[goodpix]**2)
		    RollMapPower = np.nanmean(RollMap[goodpix]**2)

		    S_2[(dz,dy+dY,dx+dX)] = (np.nanmean(D[goodpix]**2)) / (OrigMapPower + RollMapPower)
	#           S_2[(dz,-(dy+dY+1),-(dx+dX+1))] = S_2[(dz,dy+dY,dx+dX)]    
	# 	    OPTIMIZED version: unquote the above.


	# 2D Display (layer by layer)
	plt.figure(1)
	plt.subplot(121)
	plt.imshow(S_2[0], interpolation = 'none', extent = [-dX,dX,-dY,dY], vmin=0, vmax=S_2.max(), aspect='auto')
	plt.title('S_2 for dV = 0')
	plt.xlabel('dx')
	plt.ylabel('dy')

	plt.subplot(122)
	plt.imshow(S_2[dZ], interpolation = 'none', extent = [-dX,dX,-dY,dY], vmin=0, vmax=S_2.max(), aspect='auto')
	plt.title('S_2 for dV = +'+str(dZ))
	plt.xlabel('dx')
	plt.ylabel('dy')

	plt.colorbar()
	plt.savefig('image3DcubeV2.png')


	# Goal: Create a 1D plot, for each value of "dv", of the average value of structure function (inside a thin ring
	#       at radius r) versus radius. Each "dv" is a different sheet of dx,dy.
	x = np.linspace(-dX,dX,nmax)
	y = np.linspace(-dY,dY,nmax)
	xx, yy = np.meshgrid(x,y)

	maxradius = math.floor( (dX**2 + dY**2)**0.5 )
	mult = 1                        # Increases or decreases the numbers of bins used. Most accurate results at mult=1.
	reselements = math.floor(mult*maxradius)
		                        # This is the number of "resolution elements" (i.e. the number of points
		                        #      on the struct_funct vs. radius plot) that we're dealing with.
	radiusmap = (xx**2 + yy**2)**0.5

	struct_funct = np.arange(nmax*reselements).reshape(nmax,reselements)

	for i in range (0, dZ+1):
	    struct_funct[i], edges, counts = ss.binned_statistic(
		radiusmap[radiusmap<maxradius], S_2[i][radiusmap<maxradius], statistic=np.nanmean, bins = reselements)

	plt.figure(2)
	for i in range (0, dZ+1):
	    plt.plot(np.arange(reselements)/mult,struct_funct[i],label='dv='+str(i))
	plt.title('Avg. Struct. Funct. vs. Radial "Distance" from Center of S_2 Plots')
	plt.xlabel(' "Radius" ')
	plt.ylabel('Average S_2')
	plt.ylim([0.9,1.1])
	plt.legend(loc='lower right')
	plt.savefig('plot3DcubeV2.png')
