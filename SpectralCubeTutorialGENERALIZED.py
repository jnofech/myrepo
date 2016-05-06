
# 5.04.16 - Spectral Cube Tutorial (GENERALIZED)
# UPDATE (5.06.16) - Calculation for S_2 is incorrect. Do not use.

from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np


def image_make():
	cube = SpectralCube.read("paws_norot.fits")
	subcube = cube[:,400:500,400:600]

	# 1. Extracting a RECTANGULAR subcube
	# 2. Compute a moment0 map (add up in the spectral direction using `moment0 = subcube.moment(0)` Remember to take `moment0.value` for what follows.
	# 3. Calculate the structure function for a small number of offsets $\delta x =\{-1,0,1\}$ and $\delta y = \{-1,0,1\}$.  Given a map $M(x,y)$
	# 
	# $$ S_2(\delta x, \delta y) = \mathrm{mean}([M(x,y) - M(x+\delta x, y+\delta y)]^2)$$

	moment0 = subcube.moment(0,axis=0)


	dX = 5                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	nmax = abs(2*dX)+1
	S_2 = np.zeros([nmax,nmax])

	n = moment0.shape[(0)]     # The matrix "M" referred to above has n rows.
	m = moment0.shape[(1)]     # The matrix "M" referred to above has m columns.

	for dx in range (-dX,dX+1):
	    for dy in range (-dY,dY+1):
		
		M = moment0.value         # This is the matrix "M" referred to above.
		P = np.arange(n*m).reshape(n,m) # This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		D = np.zeros([n,m])   # This will be the difference between M(r) and M(r+dr).
		
		D = M - np.roll(np.roll(M,-dy,axis=0),-dx,axis=1)
		
		goodpix = (P - np.roll(P,-dy,axis=0) == -dy*m) * (P - np.roll(P,-dx,axis=1) == -dx)
		        # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		        #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		        #        that of a pixel "dy" rows below is simply dy*m.
		        # In "goodpix", pixels that have wrapped around are treated as "False".
		goodpix = goodpix.astype('float')
		goodpix[goodpix==0] = 'nan'     # Now, we can disregard the wraparound pixels entirely.
		        
		S_2[(dy+dY,dx+dX)] = (np.nanmean(D * goodpix))**2
		
	
	plt.imshow(S_2, interpolation = 'none', extent = [-dX,dX,-dY,dY])
	plt.colorbar()
	plt.xlabel('dx')
	plt.ylabel('dy')
	plt.savefig('image.png')



