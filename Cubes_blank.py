# 7.28.16 - Detects and blanks noise in a data cube.

from spectral_cube import SpectralCube
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
import os.path

print("\nWelcome to Cubes_blank! \
        \n \nAvailable functions:   \n  clean: Returns a 2D array of the RMS values of noise in\n            the galaxy's data cube, and saves a spectral\n            cube of the noise data only.")

def clean(galaxyname='M51',nmax=10):

    """
    Parameters:
    -----------
    galaxyname : string
        Name of the galaxy ('M51' or 'M33').
    vmin,vmax (UNUSED)
        The minimum and maximum channels for which
        we assume signal will appear.
    nmax : int
        Number of iterations for the noise-removing
        loop.
        
    Returns:
    --------
    RMS : numpy array
        This is a 2D array across the selected galaxy,
        
    """
    if (galaxyname=='m51') or (galaxyname=='M51'):
        filename = 'paws_norot'
    elif (galaxyname=='m33') or (galaxyname=='M33'):
        filename = 'm33.co21_iram_CLEANED'
    else:
        print "ERROR: galaxyname must be 'M51' or 'M33'."
        return
    # ^ Pick one.

    cube = SpectralCube.read(filename+".fits")
    data = cube.filled_data[:].value   # Pulls "cube"'s information (position, radial velocity, Tmax) into an array.

    p,n,m = data.shape     # The cube has n rows, m columns, and p "height units".
    RMS = np.zeros(n*m).reshape(n,m)     # An array containing the RMS values of noise brightness temperatures.
    data0 = np.copy(data)  # Copy of "data" cube, just in case.

    # Calculates RMS map, then uses this map to remove signal. The signal-less map is used
    #    to calculate the RMS map again, and so on. Iterates "nmax" times.

    for iterations in range(0,nmax):

        for j in range(0,n):
            for i in range(0,m):
                RMS[j][i] = np.nanstd(data0[:,j,i])

        # Blank regions that aren't noise:
        for k in range(0,p):
            for j in range(0,n):
                for i in range(0,m):

                    if data[k,j,i] > 3*RMS[j,i]:
                        data0[k,j,i] = np.nan
    
    # Final cube-blanking.
    
    for k in range(0,p):
        for j in range(0,n):
            for i in range(0,m):

                if data[k,j,i] > 2*RMS[j,i]:
                    data[k,j,i] = np.nan
    
    # Saves the final blanked cube. This cube should ideally contain only NOISE data. 
    if os.path.isfile(filename+'_noise.fits') == False:
        hdu = fits.PrimaryHDU(cube_s,header=cube.header)
        hdulist = fits.HDUList([hdu])
        hdulist.writeto(filename+'_noise.fits')
    else:
        print "\n ERROR: "+filename+"_noise.fits' already exists. Please delete the file and try again."
        
    return RMS
