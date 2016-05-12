
# 5.11.16 - Structure Function with Noise Correction


from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import scipy.stats as ss
import math


def image_make():
	cube = SpectralCube.read("paws_norot.fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	subcube = cube[40:80,350:550,500:700]
	noisecube = cube[100:120,350:550,500:700]

	data0 = subcube.filled_data[:]		# Overall (signal + noise) data.
	ndata0 = noisecube.filled_data[:]  	# Noise data.
	 
	dX = 100                      	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)             	# Same as above, but for "dy". For simplicity, let it be the same as dX.
	dZ = 5
	ddX = 10                       	# "Step size", for the S_2 calculation. Make sure it's a factor of dX. Leave at 1 for full calculation.
	ddY = np.copy(ddX)            	# Same as above, but for dY.
	nmax = abs(2*dX/ddX)+1
	S2 = np.zeros([dZ+1,nmax,nmax])		# This is the 2nd-order Structure Function from signal+noise.
	S2n = np.zeros([dZ+1,nmax,nmax]) 	# This is the 2nd-order Structure Function from noise only.

	p,n,m = data0.shape     # The matrix "M" referred to above has n rows, m columns, and p "height units".
	p1,n1,m1 = ndata0.shape # The matrix "M1" (basically "M", but for noise only) has n1 rows, m1 columns, and p1 "height units".


	for delz in range (0,dZ+1):                     # "delz" defined as "dz / ddZ"
	    for delx in range (-dX/ddX,dX/ddX+1):       # "delx" defined as "dx / ddX"
		for dely in range (-dY/ddY,dY/ddY+1):   # "dely" defined as "dy / ddY"
		    dx = delx*ddX
		    dy = dely*ddY
		    dz = delz				# ddZ = 1. Cannot be changed in this program.

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

		    S2[(dz,(dy+dY)/ddY,(dx+dX)/ddX)] = (np.nanmean(D[goodpix]**2))
		    
		    #_____________________________________________________________
		    
		    M1 = ndata0.value         			# This is the matrix "M" referred to above, but for noise.
		    P1 = np.arange(p1*n1*m1).reshape(p1,n1,m1) 	# This will be used to track the shifting "pixels" of M1(r) and M1(r+dr).
		    D1 = np.zeros([p1,n1,m1])   		# This will be the difference between M1(r) and M1(r+dr).

		    RollMap1 = np.roll((np.roll(np.roll(M1,-dz,axis=0),-dy,axis=1)),-dx,axis=2)
		    D1 = M1 - RollMap1

		    goodpix1 = (P1 - np.roll(P1,-dz,axis=0) == -dz*n*m) * (P1 - np.roll(P1,-dy,axis=1) == -dy*m) * (P1 - np.roll(P1,-dx,axis=2) == -dx)
		            # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		            #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		            #        that of a pixel "dy" rows below is simply dy*m.
		            # In "goodpix", pixels that have wrapped around are treated as "False".

		    S2n[(dz,(dy+dY)/ddY,(dx+dX)/ddX)] = (np.nanmean(D1[goodpix1]**2))
       


	S_2 = S2 - S2n      # This is the structure function for the Signal ONLY.		
	       
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
	plt.savefig('map_galaxyarm_nonoise.png')


	# Goal: Create a 1D plot, for each value of dz, of the average value of structure function (inside a thin ring
	#       at radius r) versus radius. Each "dz" (or, physically, dV) is a different sheet of dx,dy.
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
	#plt.ylim([0.9,1.1])
	plt.legend(loc='lower right')
	plt.savefig('plot_galaxyarm_nonoise.png')
