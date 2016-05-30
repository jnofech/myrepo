# 5.25.16 - Removes Rotation from .fits file

from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
import math
import os.path

from astropy.wcs import WCS
from reproject import reproject_interp


print('\nWelcome to Cubes_rotremoval! \n \nAvailable function: \n  cleaning: Creates a rotation-corrected version of a given .fits file.\n \n')

def cleaning(filename_cube='paws-pdbi+30m-12co10-1as.cube.fixed', filename_rotmap='paws_rotmod'):
	"""
	Argument format:
	"(filename_cube='paws-pdbi+30m-12co10-1as.cube.fixed', filename_rotmap='paws_rotmod')".
	Be sure to leave out the ".fits" extension when inputting a file name.

	filename_cube - the filename of the .fits file containing the uncorrected spectral data that we want to work with.
	filename_rotmap - the filename of the above file's corresponding rotational velocity map.

	NOTE: Currently only supports 'paws-pdbi+30m-12co10-1as.cube.fixed' (M51) and 'm33.co21_iram' (M33), along with
		their respective velocity maps.

	"""

	cube = SpectralCube.read(filename_cube+".fits")     # This is the cube of the raw, uncorrected data.
	rot = fits.open(filename_rotmap+'.fits')[0]


	# Checks if 'reprojection' is necessary. If so, then it reprojects the rotational velocity map to match
	#    the cube's dimensions.

	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.
	data0 = data.value

	if (cube.shape[1] == rot.shape[0]) and (cube.shape[2] == rot.shape[1]):
	    # The spatial dimensions of the cube and its rotation map match. Reprojection is not necessary.
	    if filename_cube =='m33.co21_iram':       	# Checks if we're using the M33 data, which is in m/s for some reason.
		array = rot.data/1000.0          	# If so, converts M33's data from m/s to km/s.
	    else:
		array = rot.data
	else:
	    # The spatial dimensions of the cube and its rotation map DO NOT match. Reprojection is necessary,
	    # 	and the cube information must appear in a moment map for it to work.
	    moment0 = cube.moment(axis=0,order=0)
	    if os.path.isfile(filename_cube+'.mom0.fits') == False:
		moment0.write(filename_cube+'.mom0.fits')
	    else:
		print "WARNING: "+filename_cube+".mom0.fits' already exists. Please delete the file and try again."
	    array, footprint = reproject_interp(rot, moment0.header)
	    if filename_cube =='m33.co21_iram':
		array = array/1000.			# Converts M33's data from m/s to km/s.


	if filename_cube =='m33.co21_iram':                                 # Checks if we're dealing with M33 or M51.
	    spec_peak = cube.spectral_axis[data0.shape[0]/2].value / 1000.  # M33's data's velocity axis is NOT centered at
		                                                            #  zero. Will be accounted for later.
	else:
	    spec_peak = 0

	velocityres = cube.header['CDELT3']
	velocityres = velocityres / 1000.		# This is the velocity resolution of the raw data file in km/s.

	# f(x) = magnitude of each entry in 3D cube, where "x" refers to the velocity we're receiving at. There's a different
	#	f(x) for each "position" on the cube.
	#
	# For each pixel on the *601x935* velocity distribution plot (which has its own *f(x)*):
	#  1. Find its *F(s)*, multiply it by *e<sup>$-i2\pi as$</sup>*, where "a" is the rotational velocity at that pixel.
	#	Note: The "s" refers to the corresponding frequencies to the Fourier Transform of f(x). This "s"
	#	should be a vector with the same length as *F(s)*, which in turn has the same length as f(x).
	#  2. Use that e<sup>$-i2\pi as$</sup> *F(s)* to find *f(x-a)*, and populate an empty matrix with the cube's dimensions with it. 
	#	When all *f(x-a)* values are found, the result SHOULD BE a 3D matrix like"cube", but corrected for rotational velocity.




	vmax,ymax,xmax = data0.shape

	cleancube = np.zeros(vmax*ymax*xmax).reshape(vmax,ymax,xmax)    # This is the final rotation-corrected cube that we'll use.

	if os.path.isfile(filename_cube+'_CLEANED.fits') == False:
		for i in range (0,xmax):
		    for j in range (0,ymax):
			fx = data0[:,j,i]      	# This is the f(x) mentioned above.
			fx_bad = np.array(np.isnan(fx), dtype=np.float)   	# This is an array of all the values in "fx" that are NaN.
			fx_temp = np.nan_to_num(fx)
		
			Fs = np.fft.fft(fx_temp)    	# This is the F(s) of f(x), including all the NaN values which were
				                    	#     (inaccurately) turned into zeroes.
			Fs_bad = np.fft.fft(fx_bad) 	# This is the F_bad(s) of fx_bad. Will be turned into F_bad(x-a) to deal with the
				                    	#     final result.
		
			s = np.fft.fftfreq(len(fx_temp)) 	# These are the frequencies corresponding to each F(s) value.
		
			shift = (array[j,i] - spec_peak) / velocityres
			phase = np.exp(2*np.pi*s*1j * shift) 	# This "j" is not to be confused with the parameter used in the For loop.
		
			fxa = np.fft.ifft(Fs*phase) 	# This is f(x-a). We just need to turn all those near-zero values back to NaN.
			fxa_bad = np.fft.ifft(Fs_bad*phase)
			fxa[fxa_bad>0.001]=np.nan
		
			cleancube[:,j,i] = fxa      	# Populates each "position" with the corrected velocity spectrum.

		hdu = fits.PrimaryHDU(cleancube,header=cube.header)
		hdulist = fits.HDUList([hdu])
		hdulist.writeto(filename_cube+'_CLEANED.fits')
	else:
		print "\n ERROR: "+filename_cube+"_CLEANED.fits' already exists. Please delete the file and try again."
