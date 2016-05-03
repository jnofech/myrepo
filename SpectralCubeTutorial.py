
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


	dX = 1                      		# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = dX                      		# Same as above, but for "dy". For simplicity, let it be the same as dX.
	n = abs(2*dX)+1
	S_2 = np.zeros([n,n])

	for dx in range (-dX,dX+1):
	    for dy in range (-dY,dY+1):
		
		M = moment0.value         	# This is the matrix "M" referred to above.
		D = np.zeros([100,100])   	# This is the difference between M(r) and M(r+dr).
		
		for y in range (0,100):
		    for x in range (0,100):
		
		        if (x+dx)>=0 and (x+dx)<100 and (y+dy)>=0 and (y+dy)<100:
		            D[(y,x)] = M[(y,x)] - M[(y+dy,x+dx)]
		        else:
		            D[(y,x)] = np.nan
		
		S_2[(dy+dY,dx+dX)] = (np.nanmean(D))**2		

	# When this whole process finishes, we now have a matrix "S_2" containing the mean of the difference squared for each value of dx and dy.


	plt.imshow(S_2)
	plt.colorbar()
	plt.xlabel('dx')
	plt.ylabel('dy')
	plt.savefig('image.png')
