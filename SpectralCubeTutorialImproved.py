
# 5.03.16 - Spectral Cube Tutorial


from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np


def image_make():
	cube = SpectralCube.read("paws_norot.fits")
	subcube = cube[:,200:300,400:500]

	# 1. Extracting a 100x100 subcube
	# 2. Compute a moment0 map (add up in the spectral direction using `moment0 = subcube.moment(0)` Remember to take `moment0.value` for what follows.
	# 3. Calculate the structure function for a small number of offsets $\delta x =\{-1,0,1\}$ and $\delta y = \{-1,0,1\}$.  Given a map $M(x,y)$
	# 
	# $$ S_2(\delta x, \delta y) = \mathrm{mean}([M(x,y) - M(x+\delta x, y+\delta y)]^2)$$

	moment0 = subcube.moment(0,axis=0)


	dX = 1                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = dX                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	n = abs(2*dX)+1
	S_2 = np.zeros([n,n])

	for dx in range (-dX,dX+1):
	    for dy in range (-dY,dY+1):
		
		M = moment0.value         # This is the matrix "M" referred to above.
		P = np.arange(100**2).reshape(100,100) # This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		D = np.zeros([100,100])   # This will be the difference between M(r) and M(r+dr).
		
		D = M - np.roll(np.roll(M,-dy,axis=0),-dx,axis=1)
		
		goodpix = (P - np.roll(P,-dy,axis=0) == -dy*100) * (P - np.roll(P,-dx,axis=1) == -dx)
		        # Here, pixels that have wrapped around are treated as "False".
		goodpix = goodpix.astype('float')
		goodpix[goodpix==0] = 'nan'     # Now, we can disregard the wraparound pixels entirely.
		        
		S_2[(dy+dY,dx+dX)] = (np.nanmean(D * goodpix))**2
		
	
	plt.imshow(S_2, interpolation = 'none', extent = [-dX,dX,-dY,dY])
	plt.colorbar()
	plt.xlabel('dx')
	plt.ylabel('dy')
	plt.savefig('image.png')



