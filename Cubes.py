
# 5.10.16 - Working with Subcubes from "paws_norot.fits", in spatial and spectral dimensions.


from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import scipy.stats as ss
import math

print('Welcome to Cubes! \n \nAvailable functions: \n  cubegen - generates a subcube. \n  structgen - generates a structure function map from a given subcube. \n  mapgen - generates some 2D maps of the structure function, for each "dv". \n  plotgen - generates a 1D plot of averaged structure function versus radius.')

def cubegen(vmin,vmax,ymin,ymax,xmin,xmax):
	"""Generates a subcube of the specified dimensions from the .fits file.

	   Argument format: "(vmin,vmax, ymin,ymax, xmin,xmax.)"
	   ^ These are the parameters of the desired subcube.
	   WARNING: Selecting too large of a subcube will hugely increase processing time.
	   If you use a large cube, be sure to set deltadeltaX to be larger in structgen."""

	cube = SpectralCube.read('paws_norot.fits')
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	subcube = cube[vmin:vmax,ymin:ymax,xmin:xmax]
	return subcube

def structgen(subcube0, deltaX=30, deltaV=3, deltadeltaX=1, normalization=True):

	"""Generates a 3D structure function map from a given subcube.
	
	   Argument format: "(<subcube>, deltaX (default=30), deltaV (default=3)),
	   deltadeltaX (default=1), normalization (default=False)".
	   "deltaX" is the maximum value of dX and dY. "deltaV" is the maximum 
	   value of dV. "deltadeltaX" is the step size along both dx and dy for the
	   structure function calculation, and it should be a factor of dX.

	   Enabling normalization will normalize the structure function within the [0,1]
	   interval. Disabling normalization will prevent this.

	   Be aware that processing time will increase with large deltaX and deltaV 
	   values, but can decrease with larger deltadeltaX at the cost of plot
	   resolution."""

	data0 = subcube0.filled_data[:]

	dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dZ = deltaV		     	# Same as above, but for "dv". We call it "dZ" here since that may make it easier to visualize.
	ddX = deltadeltaX		# "Step size", for the S_2 calculation. Make sure it's a factor of dX. Leave at 1 for full calculation.
	ddY = np.copy(ddX)            	# Same as above, but for dY.
	nmax = abs(2*dX/ddX)+1
	S_2 = np.zeros([dZ+1,nmax,nmax])

	p,n,m = data0.shape     # The subcube has n rows, m columns, and p "height units".


	for delz in range (0,dZ+1):                     # "delz" defined as "dz / ddZ"
	    for delx in range (-dX/ddX,dX/ddX+1):       # "delx" defined as "dx / ddX"
		for dely in range (-dY/ddY,dY/ddY+1):   # "dely" defined as "dy / ddY"

		    dx = delx*ddX
		    dy = dely*ddY
		    dz = delz				# ddZ = 1. Cannot be changed in this program.
		
		    M = data0.value         		# This is the subcube information (unitless).
		    P = np.arange(p*n*m).reshape(p,n,m) # This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		    D = np.zeros([p,n,m])   		# This will be the difference between M(r) and M(r+dr).

		    RollMap = np.roll((np.roll(np.roll(M,-dz,axis=0),-dy,axis=1)),-dx,axis=2)
		    D = M - RollMap

		    goodpix = (P - np.roll(P,-dz,axis=0) == -dz*n*m) * (P - np.roll(P,-dy,axis=1) == -dy*m) * (P - np.roll(P,-dx,axis=2) == -dx)
		            # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		            #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		            #        that of a pixel "dy" rows below is simply dy*m.
			    #	    Similar case for a pixel "dz" layers below.
			    #
		            # In "goodpix", pixels that have wrapped around are treated as "False".
		    if normalization==True:
			    OrigMapPower = np.nanmean(M[goodpix]**2)
			    RollMapPower = np.nanmean(RollMap[goodpix]**2)

			    S_2[(dz,(dy+dY)/ddY,(dx+dX)/ddX)] = (np.nanmean(D[goodpix]**2)) / (OrigMapPower + RollMapPower)
		    else:
			    S_2[(dz,(dy+dY)/ddY,(dx+dX)/ddX)] = (np.nanmean(D[goodpix]**2))


	return S_2


def mapgen(S_2, deltaX=30, deltaV=3):
	"""Generates and saves several 2D colour plots of the structure function versus
	   position; one for no shift in spectral dimension and one for "maximum" shift
	   in spectral dimension.

	   Argument format: (S_2, deltaX). Plots are created using the resulting 3D
	   matrix from structgen, and the same deltaX that was used in structgen."""

	dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dZ = deltaV		     	# Same as above, but for "dv". We call it "dZ" here since that may make it easier to visualize.


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

def plotgen(S_2, deltaX=30, deltaV=3, deltadeltaX=1):
	"""Generates and saves a 1D plot of the azimuthally-averaged structure function versus 
	   radius, for each value of "dv".

	   Argument format: (S_2, deltaX, deltaV, deltadeltaX). Plots are created using the 
	   resulting 3D matrix from structgen, and the same deltaX, deltaV, deltadeltaX
	   that were used in structgen."""

	# Goal: Create a 1D plot, for each value of "dv", of the average value of structure function (inside a thin ring
	#       at radius r) versus radius. Each "dv" is a different sheet of dx,dy.

	dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dZ = deltaV		     	# Same as above, but for "dv". We call it "dZ" here since that may make it easier to visualize.
	ddX = deltadeltaX
	ddY = np.copy(ddX)
	nmax = abs(2*dX/ddX)+1

	x = np.linspace(-dX/ddX,dX/ddX,nmax)
	y = np.linspace(-dY/ddY,dY/ddY,nmax)
	xx, yy = np.meshgrid(x,y)

	maxradius = ( (dX/ddX)**2 + (dY/ddY)**2 )**0.5
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
	X = (np.arange(reselements)/mult) / ((reselements-1)/mult) * (dX**2 + dY**2)**0.5
	for i in range (0, dZ+1):
	    plt.plot(X, struct_funct[i],label='dv='+str(i))
	plt.title('Avg. Struct. Funct. vs. Radial "Distance" from Center of S_2 Plots')
	plt.xlabel(' "Radius" ')
	plt.ylabel('Average S_2')
#	plt.ylim([0.9,1.1])
	plt.legend(loc='lower right')
	plt.savefig('plot3DcubeV2.png')
