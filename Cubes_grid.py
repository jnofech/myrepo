
# 6.17.16 - Activates functions from Cubes_multi.py and Cubes.py for many procedurally-generated regions.

print('\nWelcome to Cubes_grid! \n \nAvailable functions: \n  arrayM51: Activates Cubes_multi.arrayM51 for many procedurally-\n                        generated region selections in M51. \n  drawM51: Activates Cubes_multi.draw for the above subcubes.\n  arrayM33: Activates Cubes_multi.arrayM33 for many procedurally-\n                        generated region selections in M33. \n  drawM33: Activates Cubes_multi.draw for the above subcubes.\n \nThis program makes use of Cubes_multi.py and Cubes.py.\n \n')

import Cubes_multi
import Cubes_array
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from spectral_cube import SpectralCube
import astropy.units as u
import math
import scipy.stats as ss
from tempfile import TemporaryFile
from astropy.table import Table
from decimal import Decimal
import csv

def arrayM51(vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=1, deltadeltaV=1, drawmap = False, normalization=False):
	"""Activates Cubes_multi.array for many procedurally-selected regions in
	   M51, all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step sizes". Also draws maps of all regions involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=1,
	   deltadeltaV=1, drawmap=False, normalization=False).

	   WARNING: Selecting too large of a vmax-vmin will hugely increase
	   processing time."""

	galaxyname = 'M51'
	filename = "paws_norot"

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   				# Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = cube.header['DIST']			# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.


	height = 150						# Height of each selected region. Must be an even number.
	width = np.copy(height)					# Width of each selected region. Must be an even number.

	if drawmap == True:
		# Generates and saves a map of entire galaxy, with axes in units of parsecs.
		plt.figure(0)
		yshape = data.shape[1]/2.0
		xshape = data.shape[2]/2.0
		plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, extent=[-xshape*pixelwidthPC,xshape*pixelwidthPC,-yshape*pixelwidthPC,yshape*pixelwidthPC], origin='lower')
		plt.colorbar()
		fig = plt.gcf()
		fig.set_size_inches(15,7)	# Enlarges the image so as to prevent squishing.
		plt.xlabel('Distance from Centre in x-direction (pc)')
		plt.ylabel('Distance from Centre in y-direction (pc)')
		plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_entire.png')
		plt.clf()

		# Generates and saves a map of entire galaxy WITH REGIONS, with axes in units of resolution elements (for easier corresponding to filenames).
		galaxymap = plt.figure(1)
		ax1 = galaxymap.add_subplot(111)
		yshape = data.shape[1]/2.0
		xshape = data.shape[2]/2.0
		plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, origin='lower')

		for ymax in range(height, data.shape[1], height/2):
			for xmax in range(width,data.shape[2],width/2):
				ymin = ymax-height
				xmin = xmax-height
				testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
				if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
					# Checks if there are a hugely-significant number of "NaN" values in the region.
					ax1.add_patch( patches.Rectangle((xmin, ymin), width, height, fill=False, edgecolor='red'))


		fig = plt.gcf()
		fig.set_size_inches(15,7)	# Enlarges the image so as to prevent squishing.
		plt.xlabel('Resolution Units (x-direction)')
		plt.ylabel('Resolution Units (y-direction)')
		plt.colorbar()
		plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_regions_procedural.png')
		plt.clf()

	# Runs 'Cubes_multi.array(...)' for each of the regions that we're using.
	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# Checks if there are a hugely-significant number of "NaN" values in the region.
				Cubes_multi.array(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname,normalization)

def drawM51(vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=1, deltadeltaV=1, normalization=False, S2threshold=0.7):
	"""
	Activates Cubes_multi.draw and Cubes_array.generate for all of the previously-
	generated subcube selections, with the same args as arrayM51.

	The arguments MUST match the args/kwargs used in arrayM51!

	Parameters:
	-----------
	normalization : bool
		Enables or disables using the normalized S2 map
		instead of the usual one.
		If set to False, the program will not return a table
		of the S2 threshold-crossing coordinates due to the way
		`S2threshold` is handled. See below for more information.
	S2threshold : float
		This is the threshold value of S2 along the principal axis
		for which coordinates will be returned.
		That is, if we look at Cubes_array's plot of S2 along the
		principal axis, then at the point where S2 == S2threshold,
		the coordinates of that point are found and used in this
		function.
		Note that this is intended to be a percent-of-S2_max
		threshold-- but since the table of S2 threshold-crossing
		coordinates will only be used in the analysis of a
		NORMALIZED S2 map (i.e. range is roughly [0,1]), it will
		be treated simply as a number that the normalized S2 has
		to cross.
	Everything else : (various types)
		Same variables (and selected values) as in arrayM51.

	Returns:
	-----------
	t : Table
		Table displaying the cube name and the x- and y-coordinates
		of the first three S2 minima (above a certain threshold
		radius).
		Also saves the table in .csv and .bin formats, as 
		'S2_minimal_M51_(vmin)to(vmax)(_norm).csv' and
		'S2_minimal_M51_(vmin)to(vmax)(_norm).bin'.

		The "_norm" bit is added onto the end if normalization is
		activated.
	t2 : Table
		Table displaying the cube name and the x- and y-coordinates
		of the position on the NORMALIZED S2 map at which S2 crosses
		S2threshold.
		Also saves the table in .csv and .bin formats, as 
		'S2_thres_M51_(vmin)to(vmax)_norm.csv' and
		'S2_thres_M51_(vmin)to(vmax)_norm.bin'.

		The "_norm" bit is added onto the end for clarity, but this
		table will only be generated if normalization==True.
	"""

	galaxyname = 'M51'
	filename = "paws_norot"

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   				# Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = cube.header['DIST']			# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	height = 150						# Height of each selected region. Must be an even number.
	width = np.copy(height)					# Width of each selected region. Must be an even number.

	i = 0							# Counts the number of cubes that we're dealing with.

	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# ^ Checks if there are a hugely-significant number of "NaN" values in the region.
				Cubes_multi.draw(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname,normalization)
				i = i+1
	imax = i						# This is the number of cubes involved.
	i = 0							# Resets the counter.	

	cubename = [None]*imax
	ymin_array = [None]*imax
	ymax_array = [None]*imax
	xmin_array = [None]*imax
	xmax_array = [None]*imax
	xcoord1 = [None]*imax
	xcoord2 = [None]*imax
	xcoord3 = [None]*imax
	ycoord1 = [None]*imax
	ycoord2 = [None]*imax
	ycoord3 = [None]*imax
	xthres = [None]*imax
	ythres = [None]*imax

	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:	

				theta, linearray1_min, thres_radii, radlist = Cubes_array.generate(galaxyname,vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,\
													deltadeltaX,deltadeltaV,201,S2threshold, normalization)
				# ^ 'theta' is the position angle, 'radlist' are the radius values along the principal axis, 'linearray1_min' are the S_2 local minima
				#	values on this line, and 'thres_radii' are the S_2 threshold-crossing values on this line-- corresponding to 'radlist' for convenience.
				# For linearray1_min, we want to find the three closest-to-zero-but-still-above-a-threshold-radius positions along this principal axis at which there 
				#	are minima.
			
				xpositions, ypositions = extremacoords(theta,linearray1_min,radlist)		# Returns the x- and y-coordinates of three extrema near the center of the map.
				xthres[i], ythres[i] = thresholdcoords(theta,thres_radii)			# Returns the x- and y-coordinates of the radius at which S2 crosses S2threshold.
				ymin_array[i] = ymin
				ymax_array[i] = ymax
				xmin_array[i] = xmin
				xmax_array[i] = xmax

				xcoord1[i] = xpositions[0]
				xcoord2[i] = xpositions[1]
				xcoord3[i] = xpositions[2]

				ycoord1[i] = ypositions[0]
				ycoord2[i] = ypositions[1]
				ycoord3[i] = ypositions[2]

				cubename[i] = galaxyname#+"_"+str(vmin)+"to"+str(vmax)+"_"+str(ymin)+"to"+str(ymax)+"_"+str(xmin)+"to"+str(xmax)

				i = i+1

	# "t" - Table containing the regions used and the corresponding extrema coordinates.
	t = Table([cubename,ymin_array,ymax_array,xmin_array,xmax_array,ycoord1,xcoord1,ycoord2,xcoord2,ycoord3,xcoord3],names=('Cube Name','ymin','ymax','xmin','xmax',\
																'y1','x1','y2','x2','y3','x3'), meta={'name': 'TABLE'})
	t['ymin'].unit='pixels'
	t['ymax'].unit='pixels'
	t['xmin'].unit='pixels'
	t['xmax'].unit='pixels'
	t['y1'].unit='pc'
	t['y2'].unit='pc'
	t['y3'].unit='pc'
	t['x1'].unit='pc'
	t['x2'].unit='pc'
	t['x3'].unit='pc'

	# "t2" - Table containing the regions used and the corresponding coordinates of S_2 threshold locations.
	t2 = Table([cubename,ymin_array,ymax_array,xmin_array,xmax_array,ythres,xthres],names=('Cube Name','ymin','ymax','xmin','xmax','ythres','xthres'), meta={'name': 'TABLE'})

	t2['ymin'].unit='pixels'
	t2['ymax'].unit='pixels'
	t2['xmin'].unit='pixels'
	t2['xmax'].unit='pixels'
	t2['ythres'].unit='pc'
	t2['xthres'].unit='pc'

	# Save table 't' as a list in .csv format
	# Save table 't' as an array in .bin format
	if normalization==True:
		with open('S2_minimal_M51_'+str(vmin)+'to'+str(vmax)+'_norm.csv', 'w') as csvfile:	# Saves the following into 'S2_minimal_M51_40to80_norm.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t]
		f = file('S2_minimal_M51_'+str(vmin)+'to'+str(vmax)+'_norm.bin','wb')			# Saves the following into 'S2_minimal_M51_40to80_norm.bin'.
		np.save(f,t)
		f.close()
	else:
		with open('S2_minimal_M51_'+str(vmin)+'to'+str(vmax)+'.csv', 'w') as csvfile:		# Saves the following into 'S2_minimal_M51_40to80.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t]
		f = file('S2_minimal_M51_'+str(vmin)+'to'+str(vmax)+'.bin','wb')			# Saves the following into 'S2_minimal_M51_40to80.bin'.
		np.save(f,t)
		f.close()

	# Save table 't2' as a list in .csv format
	# Save table 't2' as an array in .bin format
	if normalization==True:
		with open('S2_thres_M51_'+str(vmin)+'to'+str(vmax)+'_norm.csv', 'w') as csvfile:	# Saves the following into 'S2_thres_M51_40to80_norm.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t]
		f = file('S2_thres_M51_'+str(vmin)+'to'+str(vmax)+'_norm.bin','wb')			# Saves the following into 'S2_thres_M51_40to80_norm.bin'.
		np.save(f,t)
		f.close()
	else:
		print "NOTE: Normalization must be enabled for the S2 threshold-\n \
			crossing table to be saved."							# DOESN'T save 't2' into 'S2_thres_M51_40to80.csv'.
													# DOESN'T save 't2' into 'S2_thres_M51_40to80.bin'.


	return t2

def arrayM33(vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=1, deltadeltaV=1, drawmap=False, normalization=False):
	"""Activates Cubes_multi.array for many procedurally-selected regions in
	   M33, all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step sizes". Also draws maps of all regions involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=1,
	   deltadeltaV=1, drawmap=False, normalization=False).

	   WARNING: Selecting too large of a vmax-vmin will hugely increase
	   processing time."""

	galaxyname = 'M33'
	filename = 'm33.co21_iram_CLEANED'

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = 840000.0					# The distance to the galaxy that M33's .fits file deals with, in parsecs.
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	height = 150						# Height of each selected region. Must be an even number.
	width = np.copy(height)					# Width of each selected region. Must be an even number.
	

	if drawmap == True:
		# Generates and saves a map of entire galaxy, with axes in units of parsecs.
		plt.figure(0)
		yshape = data.shape[1]/2.0
		xshape = data.shape[2]/2.0
		plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0,vmax=1, extent=[-xshape*pixelwidthPC,xshape*pixelwidthPC,-yshape*pixelwidthPC,yshape*pixelwidthPC], origin='lower')
		plt.colorbar()
		fig = plt.gcf()
		fig.set_size_inches(7,10)	# Enlarges the image so as to prevent squishing.
		plt.xlabel('Distance from Centre in x-direction (pc)')
		plt.ylabel('Distance from Centre in y-direction (pc)')
		plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_entire.png')
		plt.clf()

		# Generates and saves a map of entire galaxy WITH REGIONS, with axes in units of resolution elements (for easier corresponding to filenames).
		galaxymap = plt.figure(1)
		ax1 = galaxymap.add_subplot(111)
		yshape = data.shape[1]/2.0
		xshape = data.shape[2]/2.0
		plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0,vmax=1, origin='lower')

		for ymax in range(height, data.shape[1], height/2):
			for xmax in range(width,data.shape[2],width/2):
				ymin = ymax-height
				xmin = xmax-height
				testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
				if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
					# Checks if there are a hugely-significant number of "NaN" values in the region.
					ax1.add_patch( patches.Rectangle((xmin, ymin), width, height, fill=False, edgecolor='red'))


		fig = plt.gcf()
		fig.set_size_inches(7,10)	# Enlarges the image so as to prevent squishing.
		plt.xlabel('Resolution Units (x-direction)')
		plt.ylabel('Resolution Units (y-direction)')
		plt.colorbar()
		plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_regions_procedural.png')
		plt.clf()

	# Runs 'Cubes_multi.array(...)' for each of the regions that we're using.
	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# Checks if there are a hugely-significant number of "NaN" values in the region.
				Cubes_multi.array(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname,normalization)

def drawM33(vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=1, deltadeltaV=1,normalization=False):
	"""
	Activates Cubes_multi.draw and Cubes_array.generate for all of the previously-
	generated subcube selections, with the same args as arrayM33.

	The arguments MUST match the args/kwargs used in arrayM33!


	Returns:
	-----------
	t : Table
		Table displaying the cube name and the x- and y-coordinates
		of the first three S2 minima (above a certain threshold
		radius).
		Also saves the table in .csv and .bin formats, as 
		'S2_minimal_(vmin)to(vmax)M33(_norm).csv' and
		'S2_minimal_(vmin)to(vmax)M33(_norm).bin'.

		The "_norm" bit is added onto the end if normalization is
		activated.
	"""

	galaxyname = 'M33'
	filename = 'm33.co21_iram_CLEANED'

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = 840000.0					# The distance to the galaxy that M33's .fits file deals with, in parsecs.
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	height = 150						# Height of each selected region. Must be an even number.
	width = np.copy(height)					# Width of each selected region. Must be an even number.

	i = 0							# Counts the number of cubes that we're dealing with.

	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# ^ Checks if there are a hugely-significant number of "NaN" values in the region.
				Cubes_multi.draw(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname,normalization)
				i = i+1
	imax = i						# This is the number of cubes involved.
	i = 0							# Resets the counter.	

	cubename = [None]*imax
	ymin_array = [None]*imax
	ymax_array = [None]*imax
	xmin_array = [None]*imax
	xmax_array = [None]*imax
	xcoord1 = [None]*imax
	xcoord2 = [None]*imax
	xcoord3 = [None]*imax
	ycoord1 = [None]*imax
	ycoord2 = [None]*imax
	ycoord3 = [None]*imax

	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:	

				theta, linearray1_min, radlist = Cubes_array.generate(galaxyname,vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,201,normalization)
				# ^ 'theta' is the position angle, 'radlist' are the radius values along the principal axis, and 'linearray1_min' are the S_2 local minima
				#	values on this line-- corresponding to 'radlist' for convenience.
				# We want to find the three closest-to-zero-but-still-above-a-threshold-radius positions along this principal axis at which there are minima.
			
				xpositions, ypositions = extremacoords(theta,linearray1_min,radlist)		# Returns the x- and y-coordinates of three extrema near the center of the map.

				ymin_array[i] = ymin
				ymax_array[i] = ymax
				xmin_array[i] = xmin
				xmax_array[i] = xmax

				xcoord1[i] = xpositions[0]
				xcoord2[i] = xpositions[1]
				xcoord3[i] = xpositions[2]

				ycoord1[i] = ypositions[0]
				ycoord2[i] = ypositions[1]
				ycoord3[i] = ypositions[2]

				cubename[i] = galaxyname#+"_"+str(vmin)+"to"+str(vmax)+"_"+str(ymin)+"to"+str(ymax)+"_"+str(xmin)+"to"+str(xmax)

				i = i+1

	t = Table([cubename,ymin_array,ymax_array,xmin_array,xmax_array,ycoord1,xcoord1,ycoord2,xcoord2,ycoord3,xcoord3], names=('Cube Name','ymin','ymax','xmin','xmax','y1','x1','y2','x2','y3','x3'), meta={'name': 'TABLE'})
	t['ymin'].unit='pixels'
	t['ymax'].unit='pixels'
	t['xmin'].unit='pixels'
	t['xmax'].unit='pixels'
	t['y1'].unit='pc'
	t['y2'].unit='pc'
	t['y3'].unit='pc'
	t['x1'].unit='pc'
	t['x2'].unit='pc'
	t['x3'].unit='pc'
			
	# Save 't' as a table in .csv format
	if normalization==True:
		with open('S2_minimal_M33_'+str(vmin)+'to'+str(vmax)+'_norm.csv', 'w') as csvfile:	# Saves the following into 'S2_minimal_M33_40to80_norm.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t]
		f = file('S2_minimal_M33_'+str(vmin)+'to'+str(vmax)+'_norm.bin','wb')			# Saves the following into 'S2_minimal_M33_40to80_norm.bin'.
		np.save(f,t)
		f.close()
	else:
		with open('S2_minimal_M33_'+str(vmin)+'to'+str(vmax)+'.csv', 'w') as csvfile:		# Saves the following into 'S2_minimal_M33_40to80.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t]
		f = file('S2_minimal_M33_'+str(vmin)+'to'+str(vmax)+'.bin','wb')			# Saves the following into 'S2_minimal_M33_40to80.bin'.
		np.save(f,t)
		f.close()		
	return t


def extremacoords(theta,linearray1_min,radlist):
	radii = radlist[~np.isnan(linearray1_min)]

	minrad = 50					# Minimum radius (pc) for extrema to be considered for the table.

	rpositions = np.zeros(3)			# We'll only consider the first three minima. If we change it
							#       here, we need to change it in the "draw" functions too.
	xpositions = [None]*3
	ypositions = [None]*3
	
	if radii[radii>minrad].size > 2:
		rpositions[0] = radii[radii>minrad][0]		# The first radius above 'minrad' at which there is a minimum.
		rpositions[1] = radii[radii>minrad][1]		# The second radius above 'minrad' at which there is a minimum.
		rpositions[2] = radii[radii>minrad][2]		# The third radius above 'minrad' at which there is a minimum.
		xpositions = rpositions*np.cos(theta)			# Converts the radius values into x- and y-coordinates.
		ypositions = rpositions*np.sin(theta)

#		xpositions[0] = round(Decimal(xpositions[0]),2)	# Converts "xpositions" to Decimal format, then rounds it.
#		xpositions[1] = round(Decimal(xpositions[1]),2)	# Converts "xpositions" to Decimal format, then rounds it.
#		xpositions[2] = round(Decimal(xpositions[2]),2)	# Converts "xpositions" to Decimal format, then rounds it.

#		ypositions[0] = round(Decimal(ypositions[0]),2)	# Converts "xpositions" to Decimal format, then rounds it.
#		ypositions[1] = round(Decimal(ypositions[1]),2)	# Converts "xpositions" to Decimal format, then rounds it.
#		ypositions[2] = round(Decimal(ypositions[2]),2)	# Converts "xpositions" to Decimal format, then rounds it.

	elif radii[radii>minrad].size == 2:
		rpositions[0] = radii[radii>minrad][0]		# The first radius above 'minrad' at which there is a minimum.
		rpositions[1] = radii[radii>minrad][1]		# The second radius above 'minrad' at which there is a minimum.
		xpositions = rpositions*np.cos(theta)			# Converts the radius values into x- and y-coordinates.
		ypositions = rpositions*np.sin(theta)

#		xpositions[0] = round(Decimal(xpositions[0]),2)	# Converts "xpositions" to Decimal format, then rounds it.
#		xpositions[1] = round(Decimal(xpositions[1]),2)	# Converts "xpositions" to Decimal format, then rounds it.

#		ypositions[0] = round(Decimal(ypositions[0]),2)	# Converts "xpositions" to Decimal format, then rounds it.
#		ypositions[1] = round(Decimal(ypositions[1]),2)	# Converts "xpositions" to Decimal format, then rounds it.

	else:
		print "ERROR: Need more minima with radii above 50pc"


	return xpositions,ypositions

def thresholdcoords(theta,thres_radii)
	if thres_radii[thres_radii>0].size>0:
		xthres = thres_radii[thres_radii>0][0]*np.cos(theta)
		ythres = thres_radii[thres_radii>0][0]*np.sin(theta)
	else:
		xthres, ythres = np.nan

	return xthres, ythres

















