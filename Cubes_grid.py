
# 6.17.16 - Activates functions from Cubes(_corr)_array.py and Cubes(_corr)_multi.py for many procedurally-generated regions.

print('\nWelcome to Cubes_grid! \n \nAvailable functions: \n  arrayM51: Activates Cubes(_corr)_multi.array for many procedurally-\n                        generated region selections in M51. \n  drawM51: Activates Cubes(_corr)_multi.draw for the above subcubes.\n  arrayM33: Activates Cubes(_corr)_multi.array for many procedurally-\n                        generated region selections in M33. \n  drawM33: Activates Cubes(_corr)_multi.draw for the above subcubes.\n \nThis program makes use of Cubes(_corr)_array.py and Cubes(_corr)_multi.py.\n	Be sure to select mode="S2" or mode="xi" for each function. \n')

import Cubes_multi
import Cubes_corr_multi
import Cubes_array
import Cubes_corr_array
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

def arrayM51(mode='S2',vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=1, deltadeltaV=1, drawmap = False, normalization=False, xi_mode=0):
	"""
	Activates Cubes(_corr)_multi.array for many procedurally-selected regions in
		M51, all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
		and "step sizes". Also draws maps of all regions involved.

	Parameters:
	-----------
	mode : string
		Selects Structure Function mode ('S2'/'S_2'), which uses code
		from Cubes_multi.
		OR
		Selects Correlation Function mode ('xi'), which uses code from
		Cubes_corr_multi.
	vmin,...,deltadeltaV : int
		Parameters used in relevant S2/xi map.
		WARNING: Selecting too large of a vmax-vmin will hugely increase
		processing time.
	drawmap : bool
		Determines whether to draw and save maps of the galaxy and the
		procedurally-selected regions.
	normalization : bool
		Enables or disables using the normalized S2 map
		instead of the usual one.
		Automatically DISABLED for 'xi'. (?)
	xi_mode : int
		For xi calculations only. 
		When "xi_mode" is 0, the program will use a cube from the default 
			.fits file and a "noise cube" from that same .fits file.
		When "xi_mode" is 1, the program will use ONLY a cube from the 
			filename+"_blank" .fits file, which is assumed to have 
			NO NOISE. 
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
		if galaxyname=="M51":
			fig.set_size_inches(15,7)	# Enlarges the image so as to prevent squishing.
		else:
			fig.set_size_inches(7,10)	# Enlarges the image for M33 (which is 'taller' than it is wide).
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
				if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
					Cubes_multi.array(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname,normalization)
				elif (mode=='xi') or (mode=='Xi'):
					Cubes_corr_multi.array(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname,xi_mode)
				else:
					print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

def drawM51(mode='S2',vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=1, deltadeltaV=1, normalization=False, S2threshold=0.7, xi_mode=0):
	"""
	Activates Cubes_multi.draw and Cubes_array.generate for all of the previously-
	generated subcube selections, with the same args as arrayM51.

	The arguments MUST match the args/kwargs used in arrayM51!

	Parameters:
	-----------
	mode : string
		Selects Structure Function mode ('S2'/'S_2'), which uses code
		from Cubes_array and Cubes_multi.
		OR
		Selects Correlation Function mode ('xi'), which uses code from
		Cubes_corr_array and Cubes_corr_multi.
	normalization : bool
		Enables or disables using the normalized S2 map
		instead of the usual one.
		If set to False, the program will not return a table
		of the S2 threshold-crossing coordinates due to the way
		`S2threshold` is handled. See below for more information.
		Automatically DISABLED for 'xi'. (?)
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
		OR:
		If 'xi' mode is enabled, this is 100% - the threshold
		percentage of xi along the principal axis for which width 
		will be	measured (e.g. S2threshold=0.7 -> width will be
		measured at 30% of maximum "xi").
	xi_mode : int
		For xi calculations only. 
		When "xi_mode" is 0, the program will use a cube from the default 
			.fits file and a "noise cube" from that same .fits file.
		When "xi_mode" is 1, the program will use ONLY a cube from the 
			filename+"_blank" .fits file, which is assumed to have 
			NO NOISE.
	Everything else : (various types)
		Same variables (and selected values) as in arrayM51.

	Returns:
	-----------
	t : Table
		Table displaying the galaxy name and the x- and y-coordinates
		of the first three S2/xi minima (above a certain threshold
		radius).
		Also saves the table in .csv and .bin formats, as 
		'S2_minimal_M51_(vmin)to(vmax)(_norm).csv' and
		'S2_minimal_M51_(vmin)to(vmax)(_norm).bin'.

		The "_norm" bit is added onto the end if normalization is
		activated.
	t2 : Table
		Table displaying the galaxy name and the x- and y-coordinates
		of the position on the NORMALIZED S2/xi map at which S2/xi
		crosses	S2threshold.
		Also saves the table in .csv and .bin formats, as 
		'S2/xi_thres_M51_(vmin)to(vmax)_norm.csv' and
		'S2/xi_thres_M51_(vmin)to(vmax)_norm.bin'.

		The "_norm" bit is added onto the end for clarity, but in 'S2'
		mode, this table will only be generated if normalization==
		True.
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
				# Counts the number of usable regions.
				i = i+1

	imax = i						# This is the number of cubes involved.
	i = 0							# Resets the counter.

	coeff_a = [None]*imax					# Intercept of the log-log plot of "xi" versus radius.
	coeff_b = [None]*imax

	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# ^ Checks if there are a hugely-significant number of "NaN" values in the region.
				if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
					Cubes_multi.draw(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname,normalization)
				elif (mode=='xi') or (mode=='Xi'):
					coeff_a[i], coeff_b[i] = Cubes_corr_multi.draw(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname,xi_mode)
				else:
					print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
					return np.nan, np.nan
				print "LOOP: "+str(ymin)+", "+str(ymax)+", "+str(xmin)+", "+str(xmax)
				i = i+1
	i = 0							# Resets the counter again.

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
				if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
					theta, linearray1_min, thres_radii, radlist = Cubes_array.generate(galaxyname,vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,\
														deltadeltaX,deltadeltaV,201,S2threshold, normalization)
				elif (mode=='xi') or (mode=='Xi'):
					theta, linearray1_min, thres_radii, radlist = Cubes_corr_array.generate(galaxyname,vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,\
														deltadeltaX,deltadeltaV,201,1.0-S2threshold)
				else:
					print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

				# ^ 'theta' is the position angle, 'radlist' are the radius values along the principal axis, 'linearray1_min' are the S_2 local minima
				#	values on this line, and 'thres_radii' are the S_2 threshold-crossing values on this line-- corresponding to 'radlist' for convenience.
				# For linearray1_min, we want to find the three closest-to-zero-but-still-above-a-threshold-radius positions along this principal axis at which there 
				#	are minima.
			
				xpositions, ypositions = extremacoords(theta,linearray1_min,radlist)		# Returns the x- and y-coordinates of three extrema near the center of the map.
				xthres[i], ythres[i] = thresholdcoords(mode,theta,thres_radii, True)	# Returns the x- and y-coordinates of the radius at which S2 crosses S2threshold.

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

	# "t3" - Table containing the two coefficients of the linear fit of log(correlation function) vs log(scale) between 50pc and 250pc scales, for EACH REGION.

	if (mode=='xi') or (mode=='Xi'):
		t3 = Table([cubename,ymin_array,ymax_array,xmin_array,xmax_array,coeff_a,coeff_b],names=('Cube Name','ymin','ymax','xmin','xmax',\
													'intercept (a)', 'slope (b)'), meta={'name': 'TABLE'})
		t3['ymin'].unit='pixels'
		t3['ymax'].unit='pixels'
		t3['xmin'].unit='pixels'
		t3['xmax'].unit='pixels'


	# Save table 't' as a list in .csv format
	# Save table 't' as an array in .bin format
	if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
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
	elif (mode=='xi') or (mode=='Xi'):
		with open('xi_minimal_M51_'+str(vmin)+'to'+str(vmax)+'.csv', 'w') as csvfile:		# Saves the following into 'xi_minimal_M51_40to80.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t]
		f = file('xi_minimal_M51_'+str(vmin)+'to'+str(vmax)+'.bin','wb')			# Saves the following into 'xi_minimal_M51_40to80.bin'.
		np.save(f,t)
		f.close()
	else:
		print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."


	# Save table 't2' as a list in .csv format
	# Save table 't2' as an array in .bin format
	if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
		if normalization==True:
			with open('S2_thres_M51_'+str(vmin)+'to'+str(vmax)+'_norm.csv', 'w') as csvfile:	# Saves the following into 'S2_thres_M51_40to80_norm.csv'.
			    writer = csv.writer(csvfile)
			    [writer.writerow(r) for r in t2]
			f = file('S2_thres_M51_'+str(vmin)+'to'+str(vmax)+'_norm.bin','wb')			# Saves the following into 'S2_thres_M51_40to80_norm.bin'.
			np.save(f,t2)
			f.close()
		else:
			print "NOTE: Normalization must be enabled for the S2 threshold-\n \
				crossing table to be saved."							# DOESN'T save 't2' into 'S2_thres_M51_40to80.csv'.
														# DOESN'T save 't2' into 'S2_thres_M51_40to80.bin'.
	elif (mode=='xi') or (mode=='Xi'):
		with open('xi_thres_M51_'+str(vmin)+'to'+str(vmax)+'.csv', 'w') as csvfile:	# Saves the following into 'xi_thres_M51_40to80.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t2]
		f = file('xi_thres_M51_'+str(vmin)+'to'+str(vmax)+'.bin','wb')			# Saves the following into 'xi_thres_M51_40to80.bin'.
		np.save(f,t2)
		f.close()
	else:
		print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."


	# Save table 't3' as a list in .csv format
	# Save table 't3' as an array in .bin format
	if (mode=='xi') or (mode=='Xi'):
		with open('xi_linear_M51_'+str(vmin)+'to'+str(vmax)+'.csv', 'w') as csvfile:	# Saves the following into 'xi_linear_M51_40to80.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t3]
		f = file('xi_linear_M51_'+str(vmin)+'to'+str(vmax)+'.bin','wb')			# Saves the following into 'xi_linear_M51_40to80.bin'.
		np.save(f,t3)
		f.close()
	else:
		t3 = np.nan

	return t,t2,t3



def arrayM33(mode='S2',vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=1, deltadeltaV=1, drawmap = False, normalization=False, xi_mode=0):
	"""
	Activates Cubes(_corr)_multi.array for many procedurally-selected regions in
		M33, all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
		and "step sizes". Also draws maps of all regions involved.

	Parameters:
	-----------
	mode : string
		Selects Structure Function mode ('S2'/'S_2'), which uses code
		from Cubes_multi.
		OR
		Selects Correlation Function mode ('xi'), which uses code from
		Cubes_corr_multi.
	vmin,...,deltadeltaV : int
		Parameters used in relevant S2/xi map.
		WARNING: Selecting too large of a vmax-vmin will hugely increase
		processing time.
	drawmap : bool
		Determines whether to draw and save maps of the galaxy and the
		procedurally-selected regions.
	normalization : bool
		Enables or disables using the normalized S2 map
		instead of the usual one.
		Automatically DISABLED for 'xi'. (?)
	xi_mode : int
		For xi calculations only. 
		When "xi_mode" is 0, the program will use a cube from the default 
			.fits file and a "noise cube" from that same .fits file.
		When "xi_mode" is 1, the program will use ONLY a cube from the 
			filename+"_blank" .fits file, which is assumed to have 
			NO NOISE. 
	"""

	galaxyname = 'M33'
	filename = "m33.co21_iram_CLEANED"

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   				# Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = cube.header['DIST']			# The distance to the galaxy that M33's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
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
		if galaxyname=="M33":
			fig.set_size_inches(15,7)	# Enlarges the image so as to prevent squishing.
		else:
			fig.set_size_inches(7,10)	# Enlarges the image for M33 (which is 'taller' than it is wide).
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
				if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
					Cubes_multi.array(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname,normalization)
				elif (mode=='xi') or (mode=='Xi'):
					Cubes_corr_multi.array(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname,xi_mode)
				else:
					print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

def drawM33(mode='S2',vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=1, deltadeltaV=1, normalization=False, S2threshold=0.7, xi_mode=0):
	"""
	Activates Cubes_multi.draw and Cubes_array.generate for all of the previously-
	generated subcube selections, with the same args as arrayM33.

	The arguments MUST match the args/kwargs used in arrayM33!

	Parameters:
	-----------
	mode : string
		Selects Structure Function mode ('S2'/'S_2'), which uses code
		from Cubes_array and Cubes_multi.
		OR
		Selects Correlation Function mode ('xi'), which uses code from
		Cubes_corr_array and Cubes_corr_multi.
	normalization : bool
		Enables or disables using the normalized S2 map
		instead of the usual one.
		If set to False, the program will not return a table
		of the S2 threshold-crossing coordinates due to the way
		`S2threshold` is handled. See below for more information.
		Automatically DISABLED for 'xi'. (?)
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
		OR:
		If 'xi' mode is enabled, this is 100% - the threshold
		percentage of xi along the principal axis for which width 
		will be	measured (e.g. S2threshold=0.7 -> width will be
		measured at 30% of maximum "xi").
	xi_mode : int
		For xi calculations only. 
		When "xi_mode" is 0, the program will use a cube from the default 
			.fits file and a "noise cube" from that same .fits file.
		When "xi_mode" is 1, the program will use ONLY a cube from the 
			filename+"_blank" .fits file, which is assumed to have 
			NO NOISE.
	Everything else : (various types)
		Same variables (and selected values) as in arrayM33.

	Returns:
	-----------
	t : Table
		Table displaying the galaxy name and the x- and y-coordinates
		of the first three S2/xi minima (above a certain threshold
		radius).
		Also saves the table in .csv and .bin formats, as 
		'S2_minimal_M33_(vmin)to(vmax)(_norm).csv' and
		'S2_minimal_M33_(vmin)to(vmax)(_norm).bin'.

		The "_norm" bit is added onto the end if normalization is
		activated.
	t2 : Table
		Table displaying the galaxy name and the x- and y-coordinates
		of the position on the NORMALIZED S2/xi map at which S2/xi
		crosses	S2threshold.
		Also saves the table in .csv and .bin formats, as 
		'S2/xi_thres_M33_(vmin)to(vmax)_norm.csv' and
		'S2/xi_thres_M33_(vmin)to(vmax)_norm.bin'.

		The "_norm" bit is added onto the end for clarity, but in 'S2'
		mode, this table will only be generated if normalization==
		True.
	"""

	galaxyname = 'M33'
	filename = "m33.co21_iram_CLEANED"

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   				# Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = cube.header['DIST']			# The distance to the galaxy that M33's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
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
				# Counts the number of usable regions.
				i = i+1

	imax = i						# This is the number of cubes involved.
	i = 0							# Resets the counter.

	coeff_a = [None]*imax					# Intercept of the log-log plot of "xi" versus radius.
	coeff_b = [None]*imax

	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# ^ Checks if there are a hugely-significant number of "NaN" values in the region.
				if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
					Cubes_multi.draw(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname,normalization)
				elif (mode=='xi') or (mode=='Xi'):
					coeff_a[i], coeff_b[i] = Cubes_corr_multi.draw(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname,xi_mode)
				else:
					print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
					return np.nan, np.nan
				print "LOOP: "+str(ymin)+", "+str(ymax)+", "+str(xmin)+", "+str(xmax)
				i = i+1
	i = 0							# Resets the counter again.

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
				if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
					theta, linearray1_min, thres_radii, radlist = Cubes_array.generate(galaxyname,vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,\
														deltadeltaX,deltadeltaV,201,S2threshold, normalization)
				elif (mode=='xi') or (mode=='Xi'):
					theta, linearray1_min, thres_radii, radlist = Cubes_corr_array.generate(galaxyname,vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,\
														deltadeltaX,deltadeltaV,201,1.0-S2threshold)
				else:
					print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."

				# ^ 'theta' is the position angle, 'radlist' are the radius values along the principal axis, 'linearray1_min' are the S_2 local minima
				#	values on this line, and 'thres_radii' are the S_2 threshold-crossing values on this line-- corresponding to 'radlist' for convenience.
				# For linearray1_min, we want to find the three closest-to-zero-but-still-above-a-threshold-radius positions along this principal axis at which there 
				#	are minima.
			
				xpositions, ypositions = extremacoords(theta,linearray1_min,radlist)		# Returns the x- and y-coordinates of three extrema near the center of the map.
				xthres[i], ythres[i] = thresholdcoords(mode,theta,thres_radii, True)	# Returns the x- and y-coordinates of the radius at which S2 crosses S2threshold.

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

	# "t3" - Table containing the two coefficients of the linear fit of log(correlation function) vs log(scale) between 50pc and 250pc scales, for EACH REGION.

	if (mode=='xi') or (mode=='Xi'):
		t3 = Table([cubename,ymin_array,ymax_array,xmin_array,xmax_array,coeff_a,coeff_b],names=('Cube Name','ymin','ymax','xmin','xmax',\
													'intercept (a)', 'slope (b)'), meta={'name': 'TABLE'})
		t3['ymin'].unit='pixels'
		t3['ymax'].unit='pixels'
		t3['xmin'].unit='pixels'
		t3['xmax'].unit='pixels'


	# Save table 't' as a list in .csv format
	# Save table 't' as an array in .bin format
	if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
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
	elif (mode=='xi') or (mode=='Xi'):
		with open('xi_minimal_M33_'+str(vmin)+'to'+str(vmax)+'.csv', 'w') as csvfile:		# Saves the following into 'xi_minimal_M33_40to80.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t]
		f = file('xi_minimal_M33_'+str(vmin)+'to'+str(vmax)+'.bin','wb')			# Saves the following into 'xi_minimal_M33_40to80.bin'.
		np.save(f,t)
		f.close()
	else:
		print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."


	# Save table 't2' as a list in .csv format
	# Save table 't2' as an array in .bin format
	if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
		if normalization==True:
			with open('S2_thres_M33_'+str(vmin)+'to'+str(vmax)+'_norm.csv', 'w') as csvfile:	# Saves the following into 'S2_thres_M33_40to80_norm.csv'.
			    writer = csv.writer(csvfile)
			    [writer.writerow(r) for r in t2]
			f = file('S2_thres_M33_'+str(vmin)+'to'+str(vmax)+'_norm.bin','wb')			# Saves the following into 'S2_thres_M33_40to80_norm.bin'.
			np.save(f,t2)
			f.close()
		else:
			print "NOTE: Normalization must be enabled for the S2 threshold-\n \
				crossing table to be saved."							# DOESN'T save 't2' into 'S2_thres_M33_40to80.csv'.
														# DOESN'T save 't2' into 'S2_thres_M33_40to80.bin'.
	elif (mode=='xi') or (mode=='Xi'):
		with open('xi_thres_M33_'+str(vmin)+'to'+str(vmax)+'.csv', 'w') as csvfile:	# Saves the following into 'xi_thres_M33_40to80.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t2]
		f = file('xi_thres_M33_'+str(vmin)+'to'+str(vmax)+'.bin','wb')			# Saves the following into 'xi_thres_M33_40to80.bin'.
		np.save(f,t2)
		f.close()
	else:
		print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."


	# Save table 't3' as a list in .csv format
	# Save table 't3' as an array in .bin format
	if (mode=='xi') or (mode=='Xi'):
		with open('xi_linear_M33_'+str(vmin)+'to'+str(vmax)+'.csv', 'w') as csvfile:	# Saves the following into 'xi_linear_M33_40to80.csv'.
		    writer = csv.writer(csvfile)
		    [writer.writerow(r) for r in t3]
		f = file('xi_linear_M33_'+str(vmin)+'to'+str(vmax)+'.bin','wb')			# Saves the following into 'xi_linear_M33_40to80.bin'.
		np.save(f,t3)
		f.close()
	else:
		t3 = np.nan

	return t,t2,t3



def extremacoords(theta,linearray1_min,radlist):
	radii = radlist[~np.isnan(linearray1_min)]

	minrad = 50					# Minimum radius (pc) for extrema to be considered for the table.

	rpositions = np.zeros(3)			# We'll only consider the first three minima. If we change it
							#       here, we need to change it in the "draw" functions too.
	xpositions = np.empty(3)
	xpositions[:] = np.nan
	ypositions = np.empty(3)
	ypositions[:] = np.nan
	
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

def thresholdcoords(mode,theta,thres_radii,normalization):
	if (mode=='s2') or (mode=='S2') or (mode=='s_2') or (mode=='S_2'):
		if normalization==True:
			if thres_radii[thres_radii>0].size > 0:
				xthres = thres_radii[thres_radii>0][0]*np.cos(theta)		# The x-coord of the first radius at which S2 is above S2threshold.
				ythres = thres_radii[thres_radii>0][0]*np.sin(theta)		# The y-coord of this first radius.
			else:
				xthres, ythres = np.nan, np.nan
		else:
			xthres, ythres = np.nan, np.nan

		return xthres, ythres
	elif (mode=='xi') or (mode=='Xi'):
		if thres_radii[thres_radii>0].size > 0:
			xthres = thres_radii[thres_radii>0].max()*np.cos(theta)		# The x-coord of the highest radius at which xi/xi.max() is above 1.0-S2threshold.
			ythres = thres_radii[thres_radii>0].max()*np.sin(theta)		# The y-coord of this highest radius.
		else:
			xthres, ythres = np.nan, np.nan
		return xthres, ythres
	else:
		print "ERROR: 'mode' must be 'S2'/'S_2' or 'xi'."
		return np.nan, np.nan




















