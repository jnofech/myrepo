
# 7.07.16 - Makes a plot of S_2 minima distance from the region's centre versus the region's distance from the galactic centre. For each of the procedurally-generated grid regions.

print("\nWelcome to Cubes_grid_minima! \n \nAvailable functions: \n  plotM51: Plots extrema coordinates/distance against the location/radial \n		distance, OR plots the normalized S2 map width against\n		the radial distance from the galactic center. See\n		?Cubes_grid_minima.plotM51 for more information.\n  plotM33: Plots extrema coordinates/distance against the location/radial \n		distance, OR plots the normalized S2 map width against\n		the radial distance from the galactic center. See\n		?Cubes_grid_minima.plotM33 for more information.\n \nThis program reads the extrema positions from .bin files, which are required for\n	the functions to work. \n")

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import math
import scipy.stats as ss
from astropy.table import Table
import csv


def plotM51(vmin=40,vmax=80,deltaX=30,deltaV=3,deltadeltaX=1,deltadeltaV=1,mode=0):
    """
    mode=0 (DEFAULT) : "extrema distance mode"
        Takes the table of S_2 minima (and their coordinates) of the NORMALIZED
        S_2 surface map whose name matches the given parameters, and then 
        plots the minima distances (i.e. the distance between the minima and
        the center of the S_2 map, in parsecs) as a function of the distance 
        of each region to the centre of the galaxy (in pixels).
        
    mode=1 (UNAVAILABLE): "extrema coordinates mode"
        Takes the table of S_2 minima (and their coordinates) of the NORMALIZED
        S_2 surface map; and displays a 2D Tmax map of the galaxy, marking the 
        positions of the various regions on the map and indicating their
        respective S2 minima distances with differently-sized markers.
        
    mode=2 : "S_2 threshold mode"
        Takes the table of S_2 threshold-crossing distances (along the map's
        principal axis) of the NORMALIZED S_2 map. This threshold is set by
        `Cubes_grid.py` before the table is formed.
        Then, it plots the threshold distances (i.e. the distance between the
        S_2 threshold and the center of the S_2 map, in parsecs) as a function
        of the distance of each region to the centre of the galaxy (in pixels).

    Parameters:
    -----------
    vmin,...,deltadeltaV : int
        Parameters used in relevant S_2 map.
    mode : int (0, 1, or 2)
        Changes the mode of the function.
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
        
        f = file('S2_thres_M51_'+str(vmin)+'to'+str(vmax)+'_norm.bin','rb')
        table2 = np.load(f)
        f.close()
    else:
        table1 = np.nan
        table2 = np.nan
    
    if mode==0:
        print "Extrema Distance Mode activated"
        minima1 = np.zeros(table.size)      # Minima distance, in parsecs.
        minima2 = np.zeros(table.size)
        minima3 = np.zeros(table.size)

        rdist = np.zeros(table.size)        # Distance from the centre of galaxy to each region, in pixels.

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
        plt.ylabel("Minima Distance (pc)")
        plt.xlabel("Distance from M51's Centre (Pixels)")

        
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if normalization==True:
            plt.savefig('S2_miniplot_M51_'+str(vmin)+'to'+str(vmax)+'_norm.png')
        else:
            plt.savefig('S2_miniplot_M51_'+str(vmin)+'to'+str(vmax)+'.png')
	plt.clf()

    
    elif mode==1:
        print "Extrema Coordinates Mode activated"
        print "ERROR - Extrema Coordinates Mode unavailable."
        
        
    elif mode==2:
        print "S_2 Threshold Mode activated"
        thres_radii = np.zeros(table2.size)    # S2 Threshold-crossing distance, in parsecs.
        
        rdist = np.zeros(table2.size)        # Distance from the centre of galaxy to each region, in pixels.
        width = np.zeros(table2.size)
        
        for i in range(0,table.size):
            ymin,ymax = table2[i][1], table2[i][2]
            xmin,xmax = table2[i][3], table2[i][4]

            ythres,xthres = table2[i][5], table2[i][6]

            rdist[i] = np.sqrt( ((ymax+ymin)/2 - ycen)**2 + ((xmax+xmin)/2 - xcen)**2 )
            width[i] = np.sqrt( ythres**2 + xthres**2 )     # "Width" of structure function, i.e. radius at
                                                            #    which S2 exceeds S2threshold.
        plt.plot(rdist,width,'r.')
        plt.title("S2 Map Width vs. Region Distance from Galactic Centre")
        plt.ylabel("S2 Map Width (pc)")
        plt.xlabel("Distance from M51's Centre (Pixels)")

        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if normalization==True:
            plt.savefig('S2_thresplot_M51_'+str(vmin)+'to'+str(vmax)+'_norm.png')
        else:
            print "ERROR: Something went wrong-- normalization should be True."
	plt.clf()

    else:
        print "ERROR: Select mode=0 (Extrema Distance Mode),\
        \n              mode=1 (Extrema Coordinates Mode),\
        \n              or mode=2 (S_2 Threshold Mode)."



def plotM33(vmin=40,vmax=80,deltaX=30,deltaV=3,deltadeltaX=1,deltadeltaV=1,mode=0):
    """
    mode=0 (DEFAULT) : "extrema distance mode"
        Takes the table of S_2 minima (and their coordinates) of the NORMALIZED
        S_2 surface map whose name matches the given parameters, and then 
        plots the minima distances (i.e. the distance between the minima and
        the center of the S_2 map, in parsecs) as a function of the distance 
        of each region to the centre of the galaxy (in pixels).
        
    mode=1 (UNAVAILABLE): "extrema coordinates mode"
        Takes the table of S_2 minima (and their coordinates) of the NORMALIZED
        S_2 surface map; and displays a 2D Tmax map of the galaxy, marking the 
        positions of the various regions on the map and indicating their
        respective S2 minima distances with differently-sized markers.
        
    mode=2 : "S_2 threshold mode"
        Takes the table of S_2 threshold-crossing distances (along the map's
        principal axis) of the NORMALIZED S_2 map. This threshold is set by
        `Cubes_grid.py` before the table is formed.
        Then, it plots the threshold distances (i.e. the distance between the
        S_2 threshold and the center of the S_2 map, in parsecs) as a function
        of the distance of each region to the centre of the galaxy (in pixels).

    Parameters:
    -----------
    vmin,...,deltadeltaV : int
        Parameters used in relevant S_2 map.
    mode : int (0, 1, or 2)
        Changes the mode of the function.
    """
    
    #filename = 'paws_norot'
    filename = 'm33.co21_iram_CLEANED'

    # ^ Pick one.

    cube = SpectralCube.read(filename+".fits")
    data = cube.filled_data[:]   # Pulls "cube"'s information (position, radial velocity, Tmax) into an array.

    ycen = data[0].shape[0]/2.0     # This is the central y-value, in pixels.
    xcen = data[0].shape[1]/2.0     # This is the central x-value, in pixels.
                                    # Here, I assumed that the center of the galaxy is smack-dab
                                    #   in the middle of the 3D array.

    normalization=True             # Normalization is forcibly ENABLED.
    if normalization==True:
        f = file('S2_minimal_M33_'+str(vmin)+'to'+str(vmax)+'_norm.bin','rb')
        table = np.load(f)
        f.close()
        
        f = file('S2_thres_M33_'+str(vmin)+'to'+str(vmax)+'_norm.bin','rb')
        table2 = np.load(f)
        f.close()
    else:
        table1 = np.nan
        table2 = np.nan
    
    if mode==0:
        print "Extrema Distance Mode activated"
        minima1 = np.zeros(table.size)      # Minima distance, in parsecs.
        minima2 = np.zeros(table.size)
        minima3 = np.zeros(table.size)

        rdist = np.zeros(table.size)        # Distance from the centre of galaxy to each region, in pixels.

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
        plt.ylabel("Minima Distance (pc)")
        plt.xlabel("Distance from M33's Centre (Pixels)")

        
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if normalization==True:
            plt.savefig('S2_miniplot_M33_'+str(vmin)+'to'+str(vmax)+'_norm.png')
        else:
            plt.savefig('S2_miniplot_M33_'+str(vmin)+'to'+str(vmax)+'.png')
	plt.clf()

    
    elif mode==1:
        print "Extrema Coordinates Mode activated"
        print "ERROR - Extrema Coordinates Mode unavailable."
        
        
    elif mode==2:
        print "S_2 Threshold Mode activated"
        thres_radii = np.zeros(table2.size)    # S2 Threshold-crossing distance, in parsecs.
        
        rdist = np.zeros(table2.size)        # Distance from the centre of galaxy to each region, in pixels.
        width = np.zeros(table2.size)
        
        for i in range(0,table.size):
            ymin,ymax = table2[i][1], table2[i][2]
            xmin,xmax = table2[i][3], table2[i][4]

            ythres,xthres = table2[i][5], table2[i][6]

            rdist[i] = np.sqrt( ((ymax+ymin)/2 - ycen)**2 + ((xmax+xmin)/2 - xcen)**2 )
            width[i] = np.sqrt( ythres**2 + xthres**2 )     # "Width" of structure function, i.e. radius at
                                                            #    which S2 exceeds S2threshold.
        plt.plot(rdist,width,'r.')
        plt.title("S2 Map Width vs. Region Distance from Galactic Centre")
        plt.ylabel("S2 Map Width (pc)")
        plt.xlabel("Distance from M33's Centre (Pixels)")

        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if normalization==True:
            plt.savefig('S2_thresplot_M33_'+str(vmin)+'to'+str(vmax)+'_norm.png')
        else:
            print "ERROR: Something went wrong-- normalization should be True."
	plt.clf()

    else:
        print "ERROR: Select mode=0 (Extrema Distance Mode),\
        \n              mode=1 (Extrema Coordinates Mode),\
        \n              or mode=2 (S_2 Threshold Mode)."
