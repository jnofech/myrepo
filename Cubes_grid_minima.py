
# 7.07.16 - Makes a plot of S_2 minima distance from the region's centre versus the region's distance from the galactic centre. For each of the procedurally-generated grid regions.

print("\nWelcome to Cubes_grid_minima! \n \nAvailable functions: \n  plotM51: Plots the minima distance against the region's distance from the M51\n		centre, OR plots the positions of the minima on M51's Tmax map \n		with an indicator of minima distance; depending on the selected \n		mode.\n  plotM33: Plots the minima distance against the region's distance from the M33\n		centre, OR plots the positions of the minima on M33's Tmax map\n		with an indicator of minima distance; depending on the selected \n		mode.\n \nThis program reads the extrema positions from .bin files, which are required for\n	the functions to work. \n")

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import math
import scipy.stats as ss
from astropy.table import Table
import csv


def plotM51(vmin=40,vmax=80,deltaX=30,deltaV=3,deltadeltaX=1,deltadeltaV=1,mapmode=False):
    """
    Takes the table of S_2 minima (and their coordinates) of the NORMALIZED
        S_2 surface map whose name matches the given parameters, and then 
        plots the minima distances (i.e. the distance between the minima and
        the center of the S_2 map, in pixels) as a function of the distance 
        of each region to the center of the galaxy.
    In "mapmode" (CURRENTLY UNAVAILABLE), this function will instead display 
        a 2D Tmax map of the galaxy, showing the positions of the various 
        regions and indicating the S2 minima distances with differently-sized 
        pixels.

    Parameters:
    -----------
    vmin,...,deltadeltaV : int
        Parameters used in relevant S_2 map.
    mapmode : bool
        Enables or disables the aforementioned "mapmode".
    """
    
    filename = 'paws_norot'
    #filename = 'm33.co21_iram_CLEANED'

    # ^ Pick one.

    cube = SpectralCube.read(filename+".fits")
    data = cube.filled_data[:]   # Pulls "cube"'s information (position, radial velocity, Tmax) into an array.

    ycen = data[0].shape[0]/2.0     # This is the central y-value, in pixels.
    xcen = data[0].shape[1]/2.0     # This is the central x-value, in pixels.
                                    # Here, I assumed that the center of the galaxy is smack-dab
                                    #   in the middle of the 3D array.

    normalization=True             # Normalization is forcibly ENABLED.
    if normalization==True:
        f = file('S2_minimal_M51_'+str(vmin)+'to'+str(vmax)+'_norm.bin','rb')
        table = np.load(f)
        f.close()
    else:
        f = file('S2_minimal_M51_'+str(vmin)+'to'+str(vmax)+'.bin','rb')
        table = np.load(f)
        f.close()
    
    if mapmode==False:
        minima1 = np.zeros(table.size)      # Minima distance, in pixels.
        minima2 = np.zeros(table.size)
        minima3 = np.zeros(table.size)

        rdist = np.zeros(table.size)        # Distance from the centre of galaxy to each region, in pixels.
                                            # The pixels here and in "minima_" all represent the same distance.

        for i in range(0,table.size):
            ymin,ymax = table[i][1], table[i][2]
            xmin,xmax = table[i][3], table[i][4]

            y1,x1 = table[i][5], table[i][6]
            y2,x2 = table[i][7], table[i][8]
            y3,x3 = table[i][9], table[i][10]     

            rdist[i] = np.sqrt( ((ymax+ymin)/2 - ycen)**2 + ((xmax+xmin)/2 - xcen)**2 )
            minima1[i] = np.sqrt( y1**2 + x1**2 )
            minima2[i] = np.sqrt( y2**2 + x2**2 )
            minima3[i] = np.sqrt( y3**2 + x3**2 )

        plt.plot(rdist,minima1,'r.',rdist,minima2,'b.',rdist,minima3,'g.')
        plt.title("Extrema Distance vs. Region Distance from Galactic Centre")
        plt.xlabel("Distance from M51's Centre (Pixels)")
        plt.ylabel("Minima Distance (Pixels)")
        
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if normalization==True:
            plt.savefig('S2_miniplot_M51_'+str(vmin)+'to'+str(vmax)+'_norm.png')
        else:
            plt.savefig('S2_miniplot_M51_'+str(vmin)+'to'+str(vmax)+'.png')

        return table

    else:
        print "mapmode is on!"
