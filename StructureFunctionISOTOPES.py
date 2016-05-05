
# 5.04.16 - Spectral Cube Tutorial (GENERALIZED)


from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np


def image_make():
	cube12 = SpectralCube.read("paws-30m-12co10-23as-cube.fits")
	cube13 = SpectralCube.read("paws-30m-13co10-23as-cube.fits")
	subcube12 = cube12[:,50:120,40:60]
	subcube13 = cube13[:,50:120,40:60]

	# 1. Extracting a RECTANGULAR subcube
	# 2. Compute a moment0 map (add up in the spectral direction using `moment0 = subcube.moment(0)` Remember to take `moment0.value` for what follows.
	# 3. Calculate the structure function for a small number of offsets $\delta x =\{-1,0,1\}$ and $\delta y = \{-1,0,1\}$.  Given a map $M(x,y)$
	# 
	# $$ S_2(\delta x, \delta y) = \mathrm{mean}([M(x,y) - M(x+\delta x, y+\delta y)]^2)$$

	moment012 = subcube12.moment(0,axis=0)
	moment013 = subcube13.moment(0,axis=0)


	dX = 10                      # This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                      # Same as above, but for "dy". For simplicity, let it be the same as dX.
	nmax = abs(2*dX)+1
	S_2_12 = np.zeros([nmax,nmax])
	S_2_13 = np.zeros([nmax,nmax])

	n12,m12 = moment012.shape     # The matrix "M" referred to above has n rows and m columns.

	n13,m13 = moment013.shape


	for dx in range (-dX,dX+1):
	    for dy in range (-dY,dY+1):
		
		M12 = moment012.value         # This is the matrix "M" referred to above.
		P12 = np.arange(n12*m12).reshape(n12,m12) # This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		D12 = np.zeros([n12,m12])   # This will be the difference between M(r) and M(r+dr).
		
		RollMap12 = np.roll(np.roll(M12,-dy,axis=0),-dx,axis=1)
		D12 = M12 - RollMap12
		
		goodpix12 = (P12 - np.roll(P12,-dy,axis=0) == -dy*m12) * (P12 - np.roll(P12,-dx,axis=1) == -dx)
		        # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		        #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		        #        that of a pixel "dy" rows below is simply dy*m.
		        # In "goodpix", pixels that have wrapped around are treated as "False".
		        
		OrigMapPower12 = (np.nanmean(M12[goodpix12])**2)
		RollMapPower12 = (np.nanmean(RollMap12[goodpix12])**2)
		
		S_2_12[(dy+dY,dx+dX)] = (np.nanmean(D12[goodpix12])**2) / (OrigMapPower12 + RollMapPower12)
		
		
		M13 = moment013.value         # This is the matrix "M" referred to above.
		P13 = np.arange(n13*m13).reshape(n13,m13) # This will be used to track the shifting "pixels" of M(r) and M(r+dr).
		D13 = np.zeros([n13,m13])   # This will be the difference between M(r) and M(r+dr).
		
		RollMap13 = np.roll(np.roll(M13,-dy,axis=0),-dx,axis=1)
		D13 = M13 - RollMap13
		
		goodpix13 = (P13 - np.roll(P13,-dy,axis=0) == -dy*m13) * (P13 - np.roll(P13,-dx,axis=1) == -dx)
		        # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
		        #        pixel above or below it by exactly m. So, the difference between a pixel's value and
		        #        that of a pixel "dy" rows below is simply dy*m.
		        # In "goodpix", pixels that have wrapped around are treated as "False".
		        
		OrigMapPower13 = (np.nanmean(M13[goodpix13])**2)
		RollMapPower13 = (np.nanmean(RollMap13[goodpix13])**2)
		
		S_2_13[(dy+dY,dx+dX)] = (np.nanmean(D13[goodpix13])**2) / (OrigMapPower13 + RollMapPower13)
		
	
	plt.figure(1)
	plt.figure(1)
	plt.subplot(121)
	plt.imshow(S_2_12, interpolation = 'none', extent = [-dX,dX,-dY,dY], vmin=0, vmax=S_2_13.max(), aspect='auto')
	plt.title('S_2 for 12CO')
	plt.xlabel('dx')
	plt.ylabel('dy')

	plt.subplot(122)
	plt.imshow(S_2_13, interpolation = 'none', extent = [-dX,dX,-dY,dY], vmin=0, vmax=S_2_13.max(), aspect='auto')
	plt.title('S_2 for 13C0')
	plt.xlabel('dx')
	plt.ylabel('dy')
	plt.colorbar()
	plt.savefig('image1213.png')



