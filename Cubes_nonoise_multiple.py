
# 5.13.16 - Calculates S_2 with noise correction using functions from "Cubes.py".

print('\nWelcome to Cubes_nonoise_multiple! \n \nAvailable functions: \n  image_make: Produces a 2D map and 1D plot of the noise-corrected S_2.  \n  multiple_imagesM51: Activates image_make for several preset subcubes all at\n                        once for M51.\n  multiple_imagesM33: Activates image_make for several preset subcubes all at\n                        once for M33.\n \nThis program makes use of Cubes.py.\n \n')
import Cubes
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from spectral_cube import SpectralCube
import astropy.units as u

def image_make(vmin, vmax, ymin, ymax, xmin, xmax, deltaX = 100, deltaV = 3, deltadeltaX = 10, deltadeltaV = 1, filename="paws_norot", drawmap="False", galaxyname='M51'):
	"""
	   Generates a noise-corrected 2D map and 1D plot of S_2, from subcube of the 
	   specified dimensions; using the .fits file in "Cubes.py".

	   Argument format: "(vmin,vmax, ymin,ymax, xmin,xmax, deltaX=100, deltaV=3,
	      deltadeltaX=10, deltadeltaV=1, filename="paws_norot", drawmap="False",
	      galaxyname='M51')."
	   ^ These are the parameters of the desired subcube, along with maximum dX/dY,
	     maximum dV, "step sizes" for calculating S_2, the selected .fits file
	     name (minus the ".fits" extension), and an option to draw the subcube
	     as it appears on the galaxy T_max map (within the range specified).

	   Note that vmin and vmax will only affect the overall structure function (from
	     signal+noise), but not the noise-only structure function.


	   WARNING: Selecting too large of a subcube will hugely increase processing time.
	   If you use a large cube, be sure to set deltadeltaX/V to be larger in structgen.
	   For reliability, make sure that deltadeltaX/V is a factor of deltaX/V.

	   Be aware that processing time will increase with large deltaX and deltaV 
	   values, but can dramatically decrease with larger deltadeltaX at the cost of
	   plot resolution (or with larget deltadeltaV at the cost of the number of 1D plots)."""

	imagename = galaxyname+"_"+str(vmin)+"to"+str(vmax)+"_"+str(ymin)+"to"+str(ymax)+"_"+str(xmin)+"to"+str(xmax)

	subcube = Cubes.cubegen(vmin,vmax,ymin,ymax,xmin,xmax,filename,drawmap, imagename)	# Will draw a map of the subcube if drawmap="True".
	noisecube = Cubes.cubegen(0,20,ymin,ymax,xmin,xmax,filename,"False", imagename)		# We don't want a map of the noise, since it's just... noise. Nothing to see there.

	S2 = Cubes.structgen(subcube,deltaX,deltaV,deltadeltaX,deltadeltaV,False)
	S2n = Cubes.structgen(noisecube,deltaX,deltaV,deltadeltaX,deltadeltaV,False)
	
	S_2 = S2 - S2n		# This is the structure function from Signal ONLY.

	Cubes.mapgen(S_2, deltaX, deltaV, deltadeltaV, imagename, filename)
	Cubes.plotgen(S_2, deltaX, deltaV, deltadeltaX, deltadeltaV, imagename, filename)




def multiple_imagesM51(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=10, deltadeltaV=1):
	"""Activates image_make with each of the .py file's subcube selections,
	   all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step sizes". Be aware that this function will forcibly draw
	   all maps and plots involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=5,
	   deltadeltaV=1).

	   WARNING: Selecting too large of a vmax-vmin will hugely increase
	   processing time."""

	galaxyname = 'M51'
	filename = "paws_norot"
	drawmap = "True"

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = cube.header['DIST']			# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.


	# ----------   NEW version:  ------------

	ymin = np.array([350,200,220,350,350,100,200])	# These are the minimum "y" values of the regions that we're dealing with.
	ymax = np.array([550,400,420,550,550,300,400])	# These are the corresponding maximum "y" values of these regions.
	xmin = np.array([500,425,260,120,250,570,360])	# These are the corresponding minimum "x" values of these regions.
	xmax = np.array([700,625,460,320,450,770,560])	# These are the corresponding maximum "x" values of these regions. (Example: The first region has ymin=350, ymax=550, xmin=500, xmax=700.)
	sets = np.ravel(ymin.shape)[0]		# This is the number of regions that we're dealing with.

	# Generates and saves a map of the galaxy, with axes in units of parsecs.
	plt.figure(0)
	yshape = data.shape[1]/2.0
	xshape = data.shape[2]/2.0
	plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, extent=[-xshape*pixelwidthPC,xshape*pixelwidthPC,-yshape*pixelwidthPC,yshape*pixelwidthPC])
	plt.colorbar()
	fig = plt.gcf()
	fig.set_size_inches(15,7)	# Enlarges the image so as to prevent squishing.
	plt.xlabel('Distance from Centre in x-direction (pc)')
	plt.ylabel('Distance from Centre in y-direction (pc)')
	plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_entire.png')
	plt.clf()

	# Generates and saves a map of the galaxy WITH REGIONS, with axes in units of resolution elements (for easier corresponding to filenames).
	galaxymap = plt.figure(1)
	ax1 = galaxymap.add_subplot(111)
	yshape = data.shape[1]/2.0
	xshape = data.shape[2]/2.0
	plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0)
	for i in range(0,sets):				# Adds red rectangles highlighting the regions that we're using.
		ax1.add_patch( patches.Rectangle((xmin[i], ymin[i]), (xmax[i]-xmin[i]), (ymax[i]-ymin[i]), fill=False, edgecolor='red'))
	fig = plt.gcf()
	fig.set_size_inches(15,7)	# Enlarges the image so as to prevent squishing.
	plt.xlabel('Resolution Units (x-direction)')
	plt.ylabel('Resolution Units (y-direction)')
	plt.colorbar()
	plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_regions.png')
	plt.clf()

	# Runs 'image_make(...)' for each of the regions that we're using. For descriptions of these regions, see the "OLD" section below.
	for i in range(0,sets):
		image_make(vmin,vmax,ymin[i],ymax[i],xmin[i],xmax[i],deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)

	# ----------   OLD version:  ------------	
#	image_make(vmin,vmax,350,550,500,700,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Chunk of one of M51's arms, away from center.
#	image_make(vmin,vmax,200,400,425,625,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Galaxy arm near center.
#	image_make(vmin,vmax,220,420,260,460,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Galaxy arm near center (other side).
#	image_make(vmin,vmax,350,550,120,320,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Chunk of galaxy arm, away from center.
#	image_make(vmin,vmax,350,550,250,450,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Chunk of galaxy arm, not near center but not quite on outskirts either.
#	image_make(vmin,vmax,100,300,570,770,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Big "empty" area.
#	image_make(vmin,vmax,200,400,360,560,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Center of galaxy. (mostly centered on the nucleus)




def multiple_imagesM33(vmin=40,vmax=80, deltaX=40, deltaV=6, deltadeltaX=10, deltadeltaV=1):
	"""Activates image_make with each of the .py file's subcube selections,
	   all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step sizes". Be aware that this function will forcibly draw
	   all maps and plots involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=5,
	   deltadeltaV=1).

	   WARNING: Selecting too large of a vmax-vmin will hugely increase
	   processing time."""

	galaxyname = 'M33'
	filename = 'm33.co21_iram_CLEANED'
	drawmap = "True"

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = 840000.0					# The distance to the galaxy that M33's .fits file deals with, in parsecs.
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	# ----------   NEW version:  ------------

	ymin = np.array([350,600,650,525,300,250])	# These are the minimum "y" values of the regions that we're dealing with.
	ymax = np.array([550,800,850,725,500,450])	# These are the corresponding maximum "y" values of these regions.
	xmin = np.array([500,100,400,288,200,550])	# These are the corresponding minimum "x" values of these regions.
	xmax = np.array([700,300,600,488,400,750])	# These are the corresponding maximum "x" values of these regions. (Example: The first region has ymin=350, ymax=550, xmin=500, xmax=700.)
	sets = np.ravel(ymin.shape)[0]		# This is the number of regions that we're dealing with.
	
	# Generates and saves a map of the galaxy, with axes in units of parsecs.
	plt.figure(0)
	yshape = data.shape[1]/2.0
	xshape = data.shape[2]/2.0
	plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, vmax=1, extent=[-xshape*pixelwidthPC,xshape*pixelwidthPC,-yshape*pixelwidthPC,yshape*pixelwidthPC])
	plt.colorbar()
	fig = plt.gcf()
	fig.set_size_inches(7, 10)	# Enlarges the image so as to prevent squishing.
	plt.xlabel('Distance from Centre in x-direction (pc)')
	plt.ylabel('Distance from Centre in y-direction (pc)')
	plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_entire.png')
	plt.clf()

	# Generates and saves a map of the galaxy WITH REGIONS, with axes in units of resolution elements (for easier corresponding to filenames).
	galaxymap = plt.figure(1)
	ax1 = galaxymap.add_subplot(111)
	yshape = data.shape[1]/2.0
	xshape = data.shape[2]/2.0
	plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, vmax=1)
	for i in range(0,sets):				# Adds red rectangles highlighting the regions that we're using.
		ax1.add_patch( patches.Rectangle((xmin[i], ymin[i]), (xmax[i]-xmin[i]), (ymax[i]-ymin[i]), fill=False, edgecolor='red'))
	fig = plt.gcf()
	fig.set_size_inches(7, 10)	# Enlarges the image so as to prevent squishing.
	plt.xlabel('Resolution Units (x-direction)')
	plt.ylabel('Resolution Units (y-direction)')
	plt.colorbar()
	plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_regions.png')
	plt.clf()

	# Runs 'image_make(...)' for each of the regions that we're using. For descriptions of these regions, see the "OLD" section below.
	for i in range(0,sets):
		image_make(vmin,vmax,ymin[i],ymax[i],xmin[i],xmax[i],deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)





	
	# ----------   OLD version:  ------------
#	image_make(vmin,vmax,350,550,500,700,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Chunk containing several arms (upper right).
#	image_make(vmin,vmax,600,800,100,300,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Chunk of arm (left).
#	image_make(vmin,vmax,650,850,400,600,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Chunk several arms (lower right).
#	image_make(vmin,vmax,525,725,288,488,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Center of galaxy. (mostly centered on the nucleus)
#	image_make(vmin,vmax,300,500,200,400,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Chunk of outer arm; mostly empty (top left).
#	image_make(vmin,vmax,250,450,550,750,deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname)		# Chunk of outer arm (top right).
