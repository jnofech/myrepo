
# 7.07.16 - Makes a plot of S_2 minima distance from the region's centre versus the region's distance from the galactic centre. For each of the procedurally-generated grid regions.

print("\nWelcome to Cubes_grid_minima! \n \nAvailable functions: \n  plotM51: Plots extrema coordinates/distance against the location/radial \n		distance, OR plots the normalized S2/xi map width against\n		the radial distance from the galactic center. See\n		?Cubes_grid_minima.plotM51 for more information.\n  plotM33: Plots extrema coordinates/distance against the location/radial \n		distance, OR plots the normalized S2/xi map width against\n		the radial distance from the galactic center. See\n		?Cubes_grid_minima.plotM33 for more information.\n \nThis program reads the extrema positions from .bin files, which are required for\n	the functions to work. \n")

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import math
import scipy.stats as ss
from astropy.table import Table
import csv


def plotM51(mode='S2',vmin=40,vmax=80,deltaX=30,deltaV=3,deltadeltaX=1,deltadeltaV=1,drawmode=0):
    """
    drawmode=0 (DEFAULT) : "extrema distance mode"
        Takes the table of S_2/xi minima (and their coordinates) of the NORMALIZED
        S_2/xi surface map whose name matches the given parameters, and then 
        plots the minima distances (i.e. the distance between the minima and
        the center of the S_2/xi map, in parsecs) as a function of the distance 
        of each region to the centre of the galaxy (in pixels).
        
    drawmode=1 : "extrema coordinates mode"
        Takes the table of S_2/xi minima (and their coordinates) of the NORMALIZED
        S_2/xi surface map; and displays a 2D Tmax map of the galaxy, marking the 
        positions of the various regions on the map and indicating their
        respective S2/xi minima distances with differently-sized markers.
        
    drawmode=2 : "S_2/xi threshold mode"
        Takes the table of S_2/xi threshold-crossing distances (along the map's
        principal axis) of the xi/NORMALIZED S_2 map. This threshold is set by
        `Cubes_grid.py` before the table is formed.
        Then, it plots the threshold distances (i.e. the distance between the
        S_2/xi threshold and the center of the S_2/xi map, in parsecs) as a
        function of the distance of each region to the centre of the galaxy (in 
        pixels).
        
    drawmode=3 : "xi slope mode"
        Takes the table of coefficients of the lines-of-best-fit for xi (and the 
        locations at which these results were found); and displays a 2D Tmax map 
        of the galaxy, marking the positions of the various regions on the map 
        and indicating their respective xi best-fit slopes with differently-sized 
        markers.
	Does not work with S2.

    Parameters:
    -----------
    vmin,...,deltadeltaV : int
        Parameters used in relevant S_2/xi map.
    drawmode : int (0, 1, or 2)
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
    if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
        if normalization==True:
            f = file('S2_minimal_M51_'+str(vmin)+'to'+str(vmax)+'_norm.bin','rb')
            table = np.load(f)
            f.close()

            f = file('S2_thres_M51_'+str(vmin)+'to'+str(vmax)+'_norm.bin','rb')
            table2 = np.load(f)
            f.close()
        else:
            table = np.nan
            table2 = np.nan
    elif (mode=='xi') or (mode=='Xi'):
            f = file('xi_minimal_M51_'+str(vmin)+'to'+str(vmax)+'.bin','rb')
            table = np.load(f)
            f.close()

            f = file('xi_thres_M51_'+str(vmin)+'to'+str(vmax)+'.bin','rb')
            table2 = np.load(f)
            f.close()
            
            f = file('xi_linear_M51_'+str(vmin)+'to'+str(vmax)+'.bin','rb')
            table3 = np.load(f)
            f.close()
            
            if deltaX==0:
                table3 = np.nan    # Clears table3, so that its replacement may be used by the same code.
                
                f = file('xi_linear_M51_'+str(vmin)+'to'+str(vmax)+'_VELOCITY.bin','rb')
                table3 = np.load(f)
                f.close()
    else:
        print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
        return

    
    if drawmode==0:
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
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            plt.title("S2 Extrema Distance vs. Region Distance from Galactic Centre")
            plt.ylabel("Minima Distance (pc)")
            plt.xlabel("Distance from M51's Centre (Pixels)")
        elif (mode=='xi') or (mode=='Xi'):
            plt.title("xi Extrema Distance vs. Region Distance from Galactic Centre")
            plt.ylabel("Minima Distance (pc)")
            plt.xlabel("Distance from M51's Centre (Pixels)")
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
        plt.yscale('log')
        plt.xscale('log')


        
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            if normalization==True:
                plt.savefig('S2_miniplot_M51_'+str(vmin)+'to'+str(vmax)+'_norm.png')
            else:
                plt.savefig('S2_miniplot_M51_'+str(vmin)+'to'+str(vmax)+'.png')
        elif (mode=='xi') or (mode=='Xi'):
                plt.savefig('xi_miniplot_M51_'+str(vmin)+'to'+str(vmax)+'.png')
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."


        plt.clf()

    
    elif drawmode==1:
        print "Extrema Coordinates Mode activated"
        
        minima1 = np.zeros(table.size)      # Minima distance, in parsecs.
        minima2 = np.zeros(table.size)
        minima3 = np.zeros(table.size)
        ymin = np.zeros(table.size)
        ymax = np.zeros(table.size)
        xmin = np.zeros(table.size)
        xmax = np.zeros(table.size)

        for i in range(0,table.size):
            ymin[i],ymax[i] = table[i][1], table[i][2]
            xmin[i],xmax[i] = table[i][3], table[i][4]

            y1,x1 = table[i][5], table[i][6]
            y2,x2 = table[i][7], table[i][8]
            y3,x3 = table[i][9], table[i][10]     

            minima1[i] = np.sqrt( y1**2 + x1**2 )
            minima2[i] = np.sqrt( y2**2 + x2**2 )
            minima3[i] = np.sqrt( y3**2 + x3**2 )
            
        maxdist = max(np.nanmax(minima1),np.nanmax(minima2),np.nanmax(minima3))    # Largest measured extrema distance.
        sizemax=1000                                                # Size of the largest dot.
        
        xcoord = (xmax+xmin)/2.0
        ycoord = (ymax+ymin)/2.0
        # NOTE: In the following, the extrema distances are proportional to the RADII of their
        #   corresponding scatterplot dots.
        size3 = sizemax*(minima3/maxdist)**2
        size2 = sizemax*(minima2/maxdist)**2
        size1 = sizemax*(minima1/maxdist)**2
        
        fig, axarr = plt.subplots(nrows=1,ncols=1)
        ax1 = axarr
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        
        ax1.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, origin='lower')
        
        ax1.scatter(xcoord,ycoord,c='white',s=size3,label='3rd minima distance')
        ax1.scatter(xcoord,ycoord,c='blue',s=size2,label='2nd minima distance')
        ax1.scatter(xcoord,ycoord,c='r',s=size1,label='1st minima distance')
        ax1.scatter(xcoord,ycoord,color='k',s=0.1)
        
        ax1.set_title('Extrema Distances over Various Regions in M51')
        ax1.set_ylabel('y-position (pixels)')
        ax1.set_xlabel('x-position (pixels)')
        ax1.legend()
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            if normalization==True:
                plt.savefig('S2_minimap_M51_'+str(vmin)+'to'+str(vmax)+'_norm.png')
            else:
                print "ERROR: Something went wrong-- normalization should be True."
        elif (mode=='xi') or (mode=='Xi'):
                plt.savefig('xi_minimap_M51_'+str(vmin)+'to'+str(vmax)+'.png')
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

        plt.clf()
    
        
    elif drawmode==2:
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            print "S_2 Threshold Mode activated"
        elif (mode=='xi') or (mode=='Xi'):
            print "xi Threshold Mode activated"
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

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
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            plt.title("S2 Map Width vs. Region Distance from Galactic Centre")
            plt.ylabel("S2 Map Width (pc)")
            plt.xlabel("Distance from M51's Centre (Pixels)")
        elif (mode=='xi') or (mode=='Xi'):
            plt.title("xi Map Width vs. Region Distance from Galactic Centre")
            plt.ylabel("xi Map Width (pc)")
            plt.xlabel("Distance from M51's Centre (Pixels)")
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
        plt.yscale('log')
        plt.xscale('log')

        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            if normalization==True:
                plt.savefig('S2_thresplot_M51_'+str(vmin)+'to'+str(vmax)+'_norm.png')
            else:
                print "ERROR: Something went wrong-- normalization should be True."
        elif (mode=='xi') or (mode=='Xi'):
            plt.savefig('xi_thresplot_M51_'+str(vmin)+'to'+str(vmax)+'.png')
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

            
        plt.clf()


    elif drawmode==3:
        print "xi Slope Mode activated"
        
        if (mode!='xi') and (mode!='Xi'):
            print "ERROR: This mode only works for the mode=='xi'."
            return

        coeff_a = np.zeros(table3.size)      # Intercept of linear fit
        coeff_b = np.zeros(table3.size)      # Slope of best fit
        ymin = np.zeros(table3.size)
        ymax = np.zeros(table3.size)
        xmin = np.zeros(table3.size)
        xmax = np.zeros(table3.size)

        for i in range(0,table3.size):
            ymin[i],ymax[i] = table3[i][1], table3[i][2]
            xmin[i],xmax[i] = table3[i][3], table3[i][4]

            coeff_a[i] = table3[i][5]
            coeff_b[i] = table3[i][6]
            
        maxslope = max(coeff_b)     # Steepest positive slope.
        minslope = min(coeff_b)     # Steepest negative slope.
        
        if np.abs(maxslope) > np.abs(minslope):
            bigcircle_sign = '+'    # If the steepest slope is POSITIVE, then the largest circle will
                                    #    represent a positive slope and will be pink in colour.
        else:
            bigcircle_sign = '-'    # If the steepest slope is NEGATIVE, then the largest circle will
                                    #    represent a negative slope and will be light blue in colour.
        steepslope = max(np.abs(maxslope),np.abs(minslope))     # Abs. value of steepest slope overall.
        
        sizemax=1000                                                # Size of the largest dot.
        
        xcoord = (xmax+xmin)/2.0
        ycoord = (ymax+ymin)/2.0
        # NOTE: In the following, the extrema distances are proportional to the RADII of their
        #   corresponding scatterplot dots.
        
        
        # Dinstinguishes negative and positive slopes.
        b_pos = np.copy(coeff_b)     # Coeff_b where all negative values are np.nan.
        b_neg = np.copy(coeff_b)     # Coeff_b where all positive values are np.nan. 
        b_pos[coeff_b<0 ] = np.nan
        b_neg[coeff_b>0] = np.nan
        
        size2 = sizemax*(steepslope/steepslope)**2
        size1_pos = sizemax*(np.abs(b_pos)/steepslope)**2
        size1_neg = sizemax*(np.abs(b_neg)/steepslope)**2
            
        fig, axarr = plt.subplots(nrows=1,ncols=1)
        ax1 = axarr
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        
        ax1.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, origin='lower')
        
        # Displaying max-slope circles:
        if bigcircle_sign=='+':
            ax1.scatter(xcoord,ycoord,c='pink',s=size2,label='Steepest xi slope (positive)')
        else:
            ax1.scatter(xcoord,ycoord,c='cyan',s=size2,label='Steepest xi slope (negative)')
        
        ax1.scatter(xcoord,ycoord,c='red',s=size1_pos/2,label='xi slope (positive)')
        ax1.scatter(xcoord,ycoord,c='blue',s=size1_neg/2,label='xi slope (negative)')
        # ^ This is for getting the legend to display properly, since the rest of the dots must be
        #        looped over due to a bug.

        
        for i in range(0,xcoord.size-1):
            ax1.scatter(xcoord[i],ycoord[i],c='red',s=size1_pos[i])
            ax1.scatter(xcoord[i],ycoord[i],c='blue',s=size1_neg[i])
        # ^ The loop shouldn't be necessary; however, there is a freak bug that causes the circles to be
        #        noticeably off-centre if "s" is an array.
        ax1.scatter(xcoord,ycoord,color='k',s=0.1)
        
        if deltaX != 0:
            ax1.set_title('xi Slopes vs. Position Shift over Various Regions in M51')
        else:
            ax1.set_title('xi Slopes vs. Radial Velocity Shift over Various Regions in M51')
        ax1.set_ylabel('y-position (pixels)')
        ax1.set_xlabel('x-position (pixels)')
        ax1.legend()
        
        if deltaX != 0:
            plt.savefig('xi_slope_M51_'+str(vmin)+'to'+str(vmax)+'.png')
        else:
            plt.savefig('xi_velocityslope_M51_'+str(vmin)+'to'+str(vmax)+'.png')

	plt.clf()

    else:
        print "ERROR: Select drawmode=0 (Extrema Distance Mode),\
        \n              drawmode=1 (Extrema Coordinates Mode),\
        \n              drawmode=2 (S_2/xi Threshold Mode),\
        \n              or drawmode=3 (xi Slope Mode)."


def plotM33(mode='S2',vmin=40,vmax=80,deltaX=30,deltaV=3,deltadeltaX=1,deltadeltaV=1,drawmode=0):
    """
    drawmode=0 (DEFAULT) : "extrema distance mode"
        Takes the table of S_2/xi minima (and their coordinates) of the NORMALIZED
        S_2/xi surface map whose name matches the given parameters, and then 
        plots the minima distances (i.e. the distance between the minima and
        the center of the S_2/xi map, in parsecs) as a function of the distance 
        of each region to the centre of the galaxy (in pixels).
        
    drawmode=1 : "extrema coordinates mode"
        Takes the table of S_2/xi minima (and their coordinates) of the NORMALIZED
        S_2/xi surface map; and displays a 2D Tmax map of the galaxy, marking the 
        positions of the various regions on the map and indicating their
        respective S2/xi minima distances with differently-sized markers.
        
    drawmode=2 : "S_2/xi threshold mode"
        Takes the table of S_2/xi threshold-crossing distances (along the map's
        principal axis) of the xi/NORMALIZED S_2 map. This threshold is set by
        `Cubes_grid.py` before the table is formed.
        Then, it plots the threshold distances (i.e. the distance between the
        S_2/xi threshold and the center of the S_2/xi map, in parsecs) as a
        function of the distance of each region to the centre of the galaxy (in 
        pixels).
        
    drawmode=3 : "xi slope mode"
        Takes the table of coefficients of the lines-of-best-fit for xi (and the 
        locations at which these results were found); and displays a 2D Tmax map 
        of the galaxy, marking the positions of the various regions on the map 
        and indicating their respective xi best-fit slopes with differently-sized 
        markers.
	Does not work with S2.

    Parameters:
    -----------
    vmin,...,deltadeltaV : int
        Parameters used in relevant S_2/xi map.
    drawmode : int (0, 1, or 2)
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
    if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
        if normalization==True:
            f = file('S2_minimal_M33_'+str(vmin)+'to'+str(vmax)+'_norm.bin','rb')
            table = np.load(f)
            f.close()

            f = file('S2_thres_M33_'+str(vmin)+'to'+str(vmax)+'_norm.bin','rb')
            table2 = np.load(f)
            f.close()
        else:
            table = np.nan
            table2 = np.nan
    elif (mode=='xi') or (mode=='Xi'):
            f = file('xi_minimal_M33_'+str(vmin)+'to'+str(vmax)+'.bin','rb')
            table = np.load(f)
            f.close()

            f = file('xi_thres_M33_'+str(vmin)+'to'+str(vmax)+'.bin','rb')
            table2 = np.load(f)
            f.close()
            
            f = file('xi_linear_M33_'+str(vmin)+'to'+str(vmax)+'.bin','rb')
            table3 = np.load(f)
            f.close()
            
            if deltaX==0:
                table3 = np.nan    # Clears table3, so that its replacement may be used by the same code.
                
                f = file('xi_linear_M33_'+str(vmin)+'to'+str(vmax)+'_VELOCITY.bin','rb')
                table3 = np.load(f)
                f.close()
    else:
        print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
        return

    
    if drawmode==0:
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
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            plt.title("S2 Extrema Distance vs. Region Distance from Galactic Centre")
            plt.ylabel("Minima Distance (pc)")
            plt.xlabel("Distance from M33's Centre (Pixels)")
        elif (mode=='xi') or (mode=='Xi'):
            plt.title("xi Extrema Distance vs. Region Distance from Galactic Centre")
            plt.ylabel("Minima Distance (pc)")
            plt.xlabel("Distance from M33's Centre (Pixels)")
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
        plt.yscale('log')
        plt.xscale('log')

        
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            if normalization==True:
                plt.savefig('S2_miniplot_M33_'+str(vmin)+'to'+str(vmax)+'_norm.png')
            else:
                plt.savefig('S2_miniplot_M33_'+str(vmin)+'to'+str(vmax)+'.png')
        elif (mode=='xi') or (mode=='Xi'):
                plt.savefig('xi_miniplot_M33_'+str(vmin)+'to'+str(vmax)+'.png')
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."


        plt.clf()

    
    elif drawmode==1:
        print "Extrema Coordinates Mode activated"
        
        minima1 = np.zeros(table.size)      # Minima distance, in parsecs.
        minima2 = np.zeros(table.size)
        minima3 = np.zeros(table.size)
        ymin = np.zeros(table.size)
        ymax = np.zeros(table.size)
        xmin = np.zeros(table.size)
        xmax = np.zeros(table.size)

        for i in range(0,table.size):
            ymin[i],ymax[i] = table[i][1], table[i][2]
            xmin[i],xmax[i] = table[i][3], table[i][4]

            y1,x1 = table[i][5], table[i][6]
            y2,x2 = table[i][7], table[i][8]
            y3,x3 = table[i][9], table[i][10]     

            minima1[i] = np.sqrt( y1**2 + x1**2 )
            minima2[i] = np.sqrt( y2**2 + x2**2 )
            minima3[i] = np.sqrt( y3**2 + x3**2 )
            
        maxdist = max(np.nanmax(minima1),np.nanmax(minima2),np.nanmax(minima3))    # Largest measured extrema distance.
        sizemax=1000                                                # Size of the largest dot for M51. Set to 600 for M33.
        
        xcoord = (xmax+xmin)/2.0
        ycoord = (ymax+ymin)/2.0
        # NOTE: In the following, the extrema distances are proportional to the RADII of their
        #   corresponding scatterplot dots.
        size3 = sizemax*(minima3/maxdist)**2
        size2 = sizemax*(minima2/maxdist)**2
        size1 = sizemax*(minima1/maxdist)**2
        
        fig, axarr = plt.subplots(nrows=1,ncols=1)
        ax1 = axarr
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        
        ax1.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, origin='lower')
        
        ax1.scatter(xcoord,ycoord,c='white',s=size3,label='3rd minima distance')
        ax1.scatter(xcoord,ycoord,c='blue',s=size2,label='2nd minima distance')
        ax1.scatter(xcoord,ycoord,c='r',s=size1,label='1st minima distance')
        ax1.scatter(xcoord,ycoord,color='k',s=0.1)
        
        ax1.set_title('Extrema Distances over Various Regions in M33')
        ax1.set_ylabel('y-position (pixels)')
        ax1.set_xlabel('x-position (pixels)')
        ax1.legend()
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            if normalization==True:
                plt.savefig('S2_minimap_M33_'+str(vmin)+'to'+str(vmax)+'_norm.png')
            else:
                print "ERROR: Something went wrong-- normalization should be True."
        elif (mode=='xi') or (mode=='Xi'):
                plt.savefig('xi_minimap_M33_'+str(vmin)+'to'+str(vmax)+'.png')
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

        plt.clf()
    
        
    elif drawmode==2:
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            print "S_2 Threshold Mode activated"
        elif (mode=='xi') or (mode=='Xi'):
            print "xi Threshold Mode activated"
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

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
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            plt.title("S2 Map Width vs. Region Distance from Galactic Centre")
            plt.ylabel("S2 Map Width (pc)")
            plt.xlabel("Distance from M33's Centre (Pixels)")
        elif (mode=='xi') or (mode=='Xi'):
            plt.title("xi Map Width vs. Region Distance from Galactic Centre")
            plt.ylabel("xi Map Width (pc)")
            plt.xlabel("Distance from M33's Centre (Pixels)")
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
        plt.yscale('log')
        plt.xscale('log')

        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
            if normalization==True:
                plt.savefig('S2_thresplot_M33_'+str(vmin)+'to'+str(vmax)+'_norm.png')
            else:
                print "ERROR: Something went wrong-- normalization should be True."
        elif (mode=='xi') or (mode=='Xi'):
            plt.savefig('xi_thresplot_M33_'+str(vmin)+'to'+str(vmax)+'.png')
        else:
            print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

            
        plt.clf()


    elif drawmode==3:
        print "xi Slope Mode activated"
        
        if (mode!='xi') and (mode!='Xi'):
            print "ERROR: This mode only works for the mode=='xi'."
            return

        coeff_a = np.zeros(table3.size)      # Intercept of linear fit
        coeff_b = np.zeros(table3.size)      # Slope of best fit
        ymin = np.zeros(table3.size)
        ymax = np.zeros(table3.size)
        xmin = np.zeros(table3.size)
        xmax = np.zeros(table3.size)

        for i in range(0,table3.size):
            ymin[i],ymax[i] = table3[i][1], table3[i][2]
            xmin[i],xmax[i] = table3[i][3], table3[i][4]

            coeff_a[i] = table3[i][5]
            coeff_b[i] = table3[i][6]
            
        maxslope = max(coeff_b)     # Steepest positive slope.
        minslope = min(coeff_b)     # Steepest negative slope.
        
        if np.abs(maxslope) > np.abs(minslope):
            bigcircle_sign = '+'    # If the steepest slope is POSITIVE, then the largest circle will
                                    #    represent a positive slope and will be pink in colour.
        else:
            bigcircle_sign = '-'    # If the steepest slope is NEGATIVE, then the largest circle will
                                    #    represent a negative slope and will be light blue in colour.
        steepslope = max(np.abs(maxslope),np.abs(minslope))     # Abs. value of steepest slope overall.
        
        sizemax=600                                                # Size of the largest dot for M33.
        
        xcoord = (xmax+xmin)/2.0
        ycoord = (ymax+ymin)/2.0
        # NOTE: In the following, the extrema distances are proportional to the RADII of their
        #   corresponding scatterplot dots.
        
        
        # Dinstinguishes negative and positive slopes.
        b_pos = np.copy(coeff_b)     # Coeff_b where all negative values are np.nan.
        b_neg = np.copy(coeff_b)     # Coeff_b where all positive values are np.nan. 
        b_pos[coeff_b<0 ] = np.nan
        b_neg[coeff_b>0] = np.nan
        
        size2 = sizemax*(steepslope/steepslope)**2
        size1_pos = sizemax*(np.abs(b_pos)/steepslope)**2
        size1_neg = sizemax*(np.abs(b_neg)/steepslope)**2
            
        fig, axarr = plt.subplots(nrows=1,ncols=1)
        ax1 = axarr
        fig = plt.gcf()
        fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.
        
        ax1.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, origin='lower')
        
        # Displaying max-slope circles:
        if bigcircle_sign=='+':
            ax1.scatter(xcoord,ycoord,c='pink',s=size2,label='Steepest xi slope (positive)')
        else:
            ax1.scatter(xcoord,ycoord,c='cyan',s=size2,label='Steepest xi slope (negative)')
        
        ax1.scatter(xcoord,ycoord,color='red',s=size1_pos/2,label='xi slope (positive)')
        ax1.scatter(xcoord,ycoord,c='blue',s=size1_neg/2,label='xi slope (negative)')
        # ^ This is for getting the legend to display properly, since the rest of the dots must be
        #        looped over due to a bug.

        
        for i in range(0,xcoord.size-1):
            ax1.scatter(xcoord[i],ycoord[i],color='red',s=size1_pos[i])
            ax1.scatter(xcoord[i],ycoord[i],c='blue',s=size1_neg[i])
        # ^ The loop shouldn't be necessary; however, there is a freak bug that causes the circles to be
        #        noticeably off-centre if "s" is an array.
        ax1.scatter(xcoord,ycoord,color='k',s=0.1)
        
        if deltaX != 0:
            ax1.set_title('xi Slopes vs. Position Shift over Various Regions in M33')
        else:
            ax1.set_title('xi Slopes vs. Radial Velocity Shift over Various Regions in M33')
        ax1.set_ylabel('y-position (pixels)')
        ax1.set_xlabel('x-position (pixels)')
        ax1.legend()
        
        if deltaX != 0:
            plt.savefig('xi_slope_M33_'+str(vmin)+'to'+str(vmax)+'.png')
        else:
            plt.savefig('xi_velocityslope_M33_'+str(vmin)+'to'+str(vmax)+'.png')

	plt.clf()

    else:
        print "ERROR: Select drawmode=0 (Extrema Distance Mode),\
        \n              drawmode=1 (Extrema Coordinates Mode),\
        \n              drawmode=2 (S_2/xi Threshold Mode),\
        \n              or drawmode=3 (xi Slope Mode)."
