
# 6.17.16 - Activates functions from Cubes_multi.py and Cubes.py for many procedurally-generated regions.

print('\nWelcome to Cubes_grid! \n \nAvailable functions: \n  arrayM51: Activates Cubes_multi.S2_arrayM51 for many procedurally-\n                        generated region selections in M51. \n  drawM51: Activates Cubes_multi.S2_draw for the above subcubes.\n  arrayM33: Activates Cubes_multi.S2_arrayM33 for many procedurally-\n                        generated region selections in M33. \n  drawM51: Activates Cubes_multi.S2_draw for the above subcubes.\n \nThis program makes use of Cubes_multi.py and Cubes.py.\n \n')

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

def arrayM51(vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=10, deltadeltaV=1, drawmap = False):
	"""Activates Cubes_multi.S2_array for many procedurally-selected regions in
	   M51, all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step sizes". Also draws maps of all regions involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=10,
	   deltadeltaV=1, drawmap=False).

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

	# Runs 'Cubes_multi.S2_array(...)' for each of the regions that we're using. For descriptions of these regions, see the "OLD" section below.
	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# Checks if there are a hugely-significant number of "NaN" values in the region.
				Cubes_multi.S2_array(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)

def drawM51(vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=10, deltadeltaV=1):
	"""Activates S2_draw for all of the previously-generated subcube selections,
	   with the same args as arrayM51.

	   Argument format: "(vmin=40,vmax=80, deltaX=30, deltaV=3, deltadeltaX=10,
	   deltadeltaV=1).

	   These MUST match the args/kwargs used in arrayM51!"""

	galaxyname = 'M51'
	filename = "paws_norot"

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   				# Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = cube.header['DIST']			# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	height = 150						# Height of each selected region. Must be an even number.
	width = np.copy(height)					# Width of each selected region. Must be an even number.

	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# Checks if there are a hugely-significant number of "NaN" values in the region.
				Cubes_multi.S2_draw(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname)
				Cubes_array.generate(galaxyname,vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,201)


def arrayM33(vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=10, deltadeltaV=1, drawmap=False):
	"""Activates Cubes_multi.S2_array for many procedurally-selected regions in
	   M33, all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step sizes". Also draws maps of all regions involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=10,
	   deltadeltaV=1, drawmap=False).

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

	# Runs 'Cubes_multi.S2_array(...)' for each of the regions that we're using. For descriptions of these regions, see the "OLD" section below.
	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# Checks if there are a hugely-significant number of "NaN" values in the region.
				Cubes_multi.S2_array(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)

def drawM33(vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=10, deltadeltaV=1):
	"""Activates S2_draw for all of the previously-generated subcube selections,
	   with the same args as arrayM33.

	   Argument format: "(vmin=40,vmax=80, deltaX=30, deltaV=6, deltadeltaX=10,
	   deltadeltaV=1).

	   These MUST match the args/kwargs used in arrayM33!"""

	galaxyname = 'M33'
	filename = 'm33.co21_iram_CLEANED'

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = 840000.0					# The distance to the galaxy that M33's .fits file deals with, in parsecs.
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	height = 150						# Height of each selected region. Must be an even number.
	width = np.copy(height)					# Width of each selected region. Must be an even number.

	for ymax in range(height, data.shape[1], height/2):
		for xmax in range(width,data.shape[2],width/2):
			ymin = ymax-height
			xmin = xmax-height
			testcube = data[vmin:vmax,ymin:ymax,xmin:xmax]
			if (np.float(np.count_nonzero(np.isnan(testcube))) / np.float(np.count_nonzero(testcube))) < 0.05:
				# Checks if there are a hugely-significant number of "NaN" values in the region.
				Cubes_multi.S2_draw(vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname)
				Cubes_array.generate(galaxyname,vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV,201)
