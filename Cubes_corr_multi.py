 
# 6.15.16 - Calculates and plots xi using functions from "Cubes_corr.py".

print('\nWelcome to Cubes_corr_multi! \n \nAvailable functions: \n  array: Saves a "xi" array. \n  draw: Generates a 2D map and 1D plot of xi.  \n  arrayM51: Activates "array" for several preset subcubes all at\n                        once for M51.\n  drawM51: Activates "draw" for the above subcubes.\n  arrayM33: Activates "array" for several preset subcubes all at\n                        once for M33.\n  drawM33: Activates "draw" for the above subcubes.\n  compare_xiarray: Saves xi arrays for M51 and M33 at dV=0. \n  compare_xidraw: Draws a 1D plot comparing the above xi arrays.\n \nThis program makes use of Cubes_corr.py.\n \n')
import Cubes_corr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from spectral_cube import SpectralCube
import astropy.units as u
import math
import scipy.stats as ss
from tempfile import TemporaryFile
 
def array(vmin, vmax, ymin, ymax, xmin, xmax, deltaX = 100, deltaV = 3, deltadeltaX = 10, deltadeltaV = 1, filename="paws_norot", drawmap=False, galaxyname='M51', xi_mode=0):
	"""
	   Generates a normalized (?) array of xi, from subcube of the 
	   specified dimensions; using the .fits file in "Cubes_corr.py".

	   Argument format: "(vmin,vmax, ymin,ymax, xmin,xmax, deltaX=100, deltaV=3,
	      deltadeltaX=10, deltadeltaV=1, filename="paws_norot", drawmap=False,
	      galaxyname='M51', xi_mode=0)."
	   ^ These are the parameters of the desired subcube, along with maximum dX/dY,
	     maximum dV, "step sizes" for calculating xi, the selected .fits file
	     name (minus the ".fits" extension), and an option to draw the subcube
	     as it appears on the galaxy T_max map (within the range specified).

	   When "xi_mode" is 0, the program will use a cube from the default .fits file and
		a "noise cube" from that same .fits file.
	   When "xi_mode" is 1, the program will use ONLY a cube from the filename+"_blank"
		.fits file, which is assumed to have NO NOISE.


	   WARNING: Selecting too large of a subcube will hugely increase processing time.
	   If you use a large cube, be sure to set deltadeltaX/V to be larger in corrgen.
	   For reliability, make sure that deltadeltaX/V is a factor of deltaX/V.

	   Be aware that processing time will increase with large deltaX and deltaV 
	   values, but can dramatically decrease with larger deltadeltaX at the cost of
	   plot resolution (or with larget deltadeltaV at the cost of the number of 1D plots)."""


	imagename = galaxyname+"_"+str(vmin)+"to"+str(vmax)+"_"+str(ymin)+"to"+str(ymax)+"_"+str(xmin)+"to"+str(xmax)
	if deltadeltaX == 1 and deltadeltaV == 1:
		tempname = 'saved_xiarray_'+imagename+'_dV_is_'+str(deltaV)+'_dX_is_'+str(deltaX)+'_MAXRES'
	else:
		tempname = 'saved_xiarray_'+imagename+'_dV_is_'+str(deltaV)+'_dX_is_'+str(deltaX)

	if xi_mode==0:
		subcube = Cubes_corr.cubegen(vmin,vmax,ymin,ymax,xmin,xmax,filename,drawmap, imagename)			# Will draw a map of the subcube if drawmap=True.
		subcube_s = Cubes_corr.cubegen(vmin,vmax,ymin,ymax,xmin,xmax,filename+"_smooth",False, imagename)	# Smoothed subcube.

		xi_o = Cubes_corr.corrgen(subcube,deltaX,deltaV,deltadeltaX,deltadeltaV)
		xi_s = Cubes_corr.corrgen(subcube_s,deltaX,deltaV,deltadeltaX,deltadeltaV)

		xi = xi_o - xi_s

		# SAVE xi into an array, with the saved filename SPECIFICALLY including the parameters used.
		f = file(tempname+".bin","wb")
		np.save(f,xi)
		f.close()
	elif xi_mode==1:
		filename = filename+"_blank"
		tempname = tempname+"_blank"
		subcube = Cubes_corr.cubegen(vmin,vmax,ymin,ymax,xmin,xmax,filename,drawmap, imagename)			# Will draw a map of the subcube if drawmap=True.

		xi_o = Cubes_corr.corrgen(subcube,deltaX,deltaV,deltadeltaX,deltadeltaV)
		xi = xi_o
	else:
		print "ERROR: Select xi_mode=0 for the default calculation, or xi_mode=1 for the blanked-cube calculation."
		return
	# SAVE xi into an array, with the saved filename SPECIFICALLY including the parameters used.
	f = file(tempname+".bin","wb")
	np.save(f,xi)
	f.close()

def draw(vmin, vmax, ymin, ymax, xmin, xmax, deltaX = 100, deltaV = 3, deltadeltaX = 10, deltadeltaV = 1, filename="paws_norot", galaxyname='M51', xi_mode=0):
	"""
	   Generates plots of xi (including a 2D plot of xi vs position, and a 1D
	   plot of xi vs radius) for each "dv" from subcube of the specified 
	   dimensions; using the saved xi arrays from "array".

	   Argument format: "(vmin,vmax, ymin,ymax, xmin,xmax, deltaX=100, deltaV=3,
	      deltadeltaX=10, deltadeltaV=1, filename="paws_norot", galaxyname='M51',
	      xi_mode=0)."
	   ^ These MUST MATCH the args/kwargs used in "array"."""


	imagename = galaxyname+"_"+str(vmin)+"to"+str(vmax)+"_"+str(ymin)+"to"+str(ymax)+"_"+str(xmin)+"to"+str(xmax)
	if deltadeltaX == 1 and deltadeltaV == 1:
		tempname = 'saved_xiarray_'+imagename+'_dV_is_'+str(deltaV)+'_dX_is_'+str(deltaX)+'_MAXRES'
	else:
		tempname = 'saved_xiarray_'+imagename+'_dV_is_'+str(deltaV)+'_dX_is_'+str(deltaX)


	if xi_mode==1:
#		filename = filename+"_blank"
		tempname = tempname+"_blank"


	# File-loading.
	f = file(tempname+".bin","rb")
	xi = np.load(f)
	f.close()

	coeff_a, coeff_b = Cubes_corr.plotgen(xi, deltaX, deltaV, deltadeltaX, deltadeltaV, imagename, filename)		# a = intercept, b = slope
	if deltaX != 0:
		Cubes_corr.mapgen(xi, deltaX, deltaV, deltadeltaV, imagename, filename)
		Cubes_corr.everythinggen(vmin, vmax, ymin, ymax, xmin, xmax, xi, deltaX, deltaV, deltadeltaX, deltadeltaV, imagename, filename)

	return coeff_a, coeff_b


def arrayM51(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=1, deltadeltaV=1, drawmap = False, xi_mode=0):
	"""Activates "array" for M51 with each of the .py file's subcube selections,
	   all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step sizes". Also draws maps of all regions involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=10,
	   deltadeltaV=1, drawmap=False, xi_mode=0).

	   When "xi_mode" is 0, the program will use a cube from the default .fits file and
		a "noise cube" from that same .fits file.
	   When "xi_mode" is 1, the program will use ONLY a cube from the filename+"_blank"
		.fits file, which is assumed to have NO NOISE.

	   WARNING: Selecting too large of a vmax-vmin will hugely increase
	   processing time."""

	galaxyname = 'M51'
	filename = "paws_norot"

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = cube.header['DIST']			# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.


	ymin = np.array([350,200,220,350,350,100,200])	# These are the minimum "y" values of the regions that we're dealing with.
	ymax = np.array([550,400,420,550,550,300,400])	# These are the corresponding maximum "y" values of these regions.
	xmin = np.array([500,425,260,120,250,570,360])	# These are the corresponding minimum "x" values of these regions.
	xmax = np.array([700,625,460,320,450,770,560])	# These are the corresponding maximum "x" values of these regions. (Example: The first region has ymin=350, ymax=550, xmin=500, xmax=700.)
	sets = np.ravel(ymin.shape)[0]		# This is the number of regions that we're dealing with.

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
		for i in range(0,sets):				# Adds red rectangles highlighting the regions that we're using.
			ax1.add_patch( patches.Rectangle((xmin[i], ymin[i]), (xmax[i]-xmin[i]), (ymax[i]-ymin[i]), fill=False, edgecolor='red'))
		fig = plt.gcf()
		fig.set_size_inches(15,7)	# Enlarges the image so as to prevent squishing.
		plt.xlabel('Resolution Units (x-direction)')
		plt.ylabel('Resolution Units (y-direction)')
		plt.colorbar()
		plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_regions.png')
		plt.clf()

	# Runs 'array(...)' for each of the regions that we're using. For descriptions of these regions, see the "OLD" section below.
	for i in range(0,sets):
		array(vmin,vmax,ymin[i],ymax[i],xmin[i],xmax[i],deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname,xi_mode)
		
def drawM51(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=1, deltadeltaV=1,xi_mode=0):
	"""Activates "draw" with each of the .py file's subcube selections,
	   with the same args as "arrayM51".

	   Argument format: "(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=10,
	   deltadeltaV=1, xi_mode=0).

	   These MUST match the args/kwargs used in "arrayM51"!"""

	galaxyname = 'M51'
	filename = "paws_norot"

	cube = SpectralCube.read(filename+".fits")

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = cube.header['DIST']			# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.


	ymin = np.array([350,200,220,350,350,100,200])	# These are the minimum "y" values of the regions that we're dealing with.
	ymax = np.array([550,400,420,550,550,300,400])	# These are the corresponding maximum "y" values of these regions.
	xmin = np.array([500,425,260,120,250,570,360])	# These are the corresponding minimum "x" values of these regions.
	xmax = np.array([700,625,460,320,450,770,560])	# These are the corresponding maximum "x" values of these regions. (Example: The first region has ymin=350, ymax=550, xmin=500, xmax=700.)
	sets = np.ravel(ymin.shape)[0]		# This is the number of regions that we're dealing with.

	for i in range(0,sets):
		draw(vmin,vmax,ymin[i],ymax[i],xmin[i],xmax[i],deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname,xi_mode)



def arrayM33(vmin=40,vmax=80, deltaX=40, deltaV=6, deltadeltaX=1, deltadeltaV=1, drawmap=False, xi_mode=0):
	"""Activates "array" for M33 with each of the .py file's subcube selections,
	   all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step sizes". Also draws maps of all regions involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=40, deltaV=6, deltadeltaX=10,
	   deltadeltaV=1, drawmap=False, xi_mode=0).

	   When "xi_mode" is 0, the program will use a cube from the default .fits file and
		a "noise cube" from that same .fits file.
	   When "xi_mode" is 1, the program will use ONLY a cube from the filename+"_blank"
		.fits file, which is assumed to have NO NOISE.

	   WARNING: Selecting too large of a vmax-vmin will hugely increase
	   processing time."""

	galaxyname = 'M33'
	filename = 'm33.co21_iram_CLEANED'

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = 840000.0					# The distance to the galaxy that M33's .fits file deals with, in parsecs.
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.



	ymin = np.array([350,600,650,525,300,250])	# These are the minimum "y" values of the regions that we're dealing with.
	ymax = np.array([550,800,850,725,500,450])	# These are the corresponding maximum "y" values of these regions.
	xmin = np.array([500,100,400,288,200,550])	# These are the corresponding minimum "x" values of these regions.
	xmax = np.array([700,300,600,488,400,750])	# These are the corresponding maximum "x" values of these regions. (Example: The first region has ymin=350, ymax=550, xmin=500, xmax=700.)
	sets = np.ravel(ymin.shape)[0]		# This is the number of regions that we're dealing with.
	

	if drawmap == True:
		# Generates and saves a map of entire galaxy, with axes in units of parsecs.
		plt.figure(0)
		yshape = data.shape[1]/2.0
		xshape = data.shape[2]/2.0
		plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, vmax=1, extent=[-xshape*pixelwidthPC,xshape*pixelwidthPC,-yshape*pixelwidthPC,yshape*pixelwidthPC], origin='lower')
		plt.colorbar()
		fig = plt.gcf()
		fig.set_size_inches(7, 10)	# Enlarges the image so as to prevent squishing.
		plt.xlabel('Distance from Centre in x-direction (pc)')
		plt.ylabel('Distance from Centre in y-direction (pc)')
		plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_entire.png')
		plt.clf()

		# Generates and saves a map of entire galaxy WITH REGIONS, with axes in units of resolution elements (for easier corresponding to filenames).
		galaxymap = plt.figure(1)
		ax1 = galaxymap.add_subplot(111)
		yshape = data.shape[1]/2.0
		xshape = data.shape[2]/2.0
		plt.imshow(np.nanmax(data[vmin:vmax].value,axis=0), vmin=0, vmax=1, origin='lower')
		for i in range(0,sets):				# Adds red rectangles highlighting the regions that we're using.
			ax1.add_patch( patches.Rectangle((xmin[i], ymin[i]), (xmax[i]-xmin[i]), (ymax[i]-ymin[i]), fill=False, edgecolor='red'))
		fig = plt.gcf()
		fig.set_size_inches(7, 10)	# Enlarges the image so as to prevent squishing.
		plt.xlabel('Resolution Units (x-direction)')
		plt.ylabel('Resolution Units (y-direction)')
		plt.colorbar()
		plt.savefig('galaxy_'+galaxyname+'_'+str(vmin)+'to'+str(vmax)+'_regions.png')
		plt.clf()

	# Runs 'array(...)' for each of the regions that we're using. For descriptions of these regions, see the "OLD" section below.
	for i in range(0,sets):
		array(vmin,vmax,ymin[i],ymax[i],xmin[i],xmax[i],deltaX,deltaV,deltadeltaX,deltadeltaV,filename,drawmap,galaxyname,xi_mode)

def drawM33(vmin=40,vmax=80, deltaX=40, deltaV=6, deltadeltaX=1, deltadeltaV=1,xi_mode=0):
	"""Activates "draw" with each of the .py file's subcube selections,
	   with the same args as "arrayM33".

	   Argument format: "(vmin=40,vmax=80, deltaX=40, deltaV=6, deltadeltaX=10,
	   deltadeltaV=1, xi_mode=0).

	   These MUST match the args/kwargs used in "arrayM33"!"""

	galaxyname = 'M33'
	filename = 'm33.co21_iram_CLEANED'

	cube = SpectralCube.read(filename+".fits")

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	distancePC = 840000.0			# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.


	ymin = np.array([350,600,650,525,300,250])	# These are the minimum "y" values of the regions that we're dealing with.
	ymax = np.array([550,800,850,725,500,450])	# These are the corresponding maximum "y" values of these regions.
	xmin = np.array([500,100,400,288,200,550])	# These are the corresponding minimum "x" values of these regions.
	xmax = np.array([700,300,600,488,400,750])	# These are the corresponding maximum "x" values of these regions. (Example: The first region has ymin=350, ymax=550, xmin=500, xmax=700.)
	sets = np.ravel(ymin.shape)[0]		# This is the number of regions that we're dealing with.

	for i in range(0,sets):
		draw(vmin,vmax,ymin[i],ymax[i],xmin[i],xmax[i],deltaX,deltaV,deltadeltaX,deltadeltaV,filename,galaxyname,xi_mode)



def compare_xiarray(vmin=40,vmax=80, deltaX=40, deltadeltaX=10):
	"""
	(NO LONGER SUPPORTED)
	For both M51 and M33, saves a numpy array containing xi as a
	   function of dx, dy, and dv (where dv is only zero). These arrays
	   can be called on in 'compare_xidraw'.
	   Argument format: "(vmin=40,vmax=80, deltaX=40, deltadeltaX=5)."

	   WARNING: Selecting too large of a vmax-vmin will hugely increase
	   processing time."""


	galaxynameM51 = 'M51'
	galaxynameM33 = 'M33'
	filenameM51 = 'paws_norot'
	filenameM33 = 'm33.co21_iram_CLEANED'
	drawmap = False
	if deltadeltaX == 1:
		savename = 'savedfile_M51andM33_comparison_MAXRES'
	else:
		savename = 'savedfile_M51andM33_comparison'

	yminM51 = 200			# 'ymin' for M51.
	ymaxM51 = 400			# etc.
	xminM51 = 360
	xmaxM51 = 560
	
	yminM33 = 525			# 'ymin' for M33.
	ymaxM33 = 725			# etc.
	xminM33 = 288
	xmaxM33 = 488


	subcubeM51 = Cubes_corr.cubegen(vmin,vmax,yminM51,ymaxM51,xminM51,xmaxM51, filenameM51, drawmap)	# Subcube for M51.


	subcubeM33 = Cubes_corr.cubegen(vmin,vmax,yminM33,ymaxM33,xminM33,xmaxM33, filenameM33, drawmap)	# Subcube for M33.


	xi_M51 = Cubes_corr.corrgen(subcubeM51,deltaX,0,deltadeltaX,1)
	xi_M33 = Cubes_corr.corrgen(subcubeM33,deltaX,0,deltadeltaX,1)



	# File-saving.
	f = file(savename+".bin","wb")
	np.save(f,xi_M51)
	np.save(f,xi_M33)
	f.close()


def compare_xidraw(vmin=40,vmax=80, deltaX=40, deltadeltaX=10):
	""" 
	(NO LONGER SUPPORTED)
	Using the saved xi arrays from 'compare_xigen', creates a plot
	of xi versus radius for M51 and M33.

	Argument format: (vmin=40, vmax=80, deltaX=40, deltadeltaX=10).
	vmin, vmax are parameters of the desired subcube.
	deltaX is the maximum value of dX that the xi arrays were found over.
	deltadeltaX is the "step size", and should be a factor of deltaX.

	All four arguments should be the same one used in compare_xigen."""

	# Goal: Create a 1D plot, for dv=0, of the average value of correlation function (inside a thin ring
	#       at radius r) versus radius. The values for M51 and M33 are on the same plot.

	if deltadeltaX == 1:
		savename = 'savedfile_M51andM33_comparison_MAXRES'
	else:
		savename = 'savedfile_M51andM33_comparison'

	# File-loading.
	f = file(savename+".bin","rb")
	xi_M51 = np.load(f)
	xi_M33 = np.load(f)
	f.close()

	dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dV = 0			     	# Same as above, but for "dv".
	ddX = deltadeltaX
	ddY = np.copy(ddX)
	ddV = 1
	nmax = abs(2*dX/ddX)+1

	filenameM51 = 'paws_norot'
	filenameM33 = 'm33.co21_iram_CLEANED'

	cubeM51 = SpectralCube.read(filenameM51+".fits")
	cubeM33 = SpectralCube.read(filenameM33+".fits")
	pixelwidthDEG_M51 = cubeM51.header['CDELT2']	# The width of each pixel, in degrees.
	pixelwidthDEG_M33 = cubeM33.header['CDELT2']	# The width of each pixel, in degrees.
	distancePC_M51 = cubeM51.header['DIST']			# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	distancePC_M33 = 840000.0				# The distance to the galaxy that M33's .fits file deals with, in parsecs. ONLY works on the CLEANED file!
	pixelwidthPC_M51 = pixelwidthDEG_M51*np.pi/180.0*distancePC_M51		# The width of each pixel, in pc.
	pixelwidthPC_M33 = pixelwidthDEG_M33*np.pi/180.0*distancePC_M33		# The width of each pixel, in pc.
	velocityresM51 = cubeM51.header['CDELT3']		# Velocity resolution in km/s.
	velocityresM33 = cubeM33.header['CDELT3'] / 1000.0	# Velocity resolution in km/s. (Header in m/s by default for M33.)

	x = np.linspace(-dX/ddX,dX/ddX,nmax)
	y = np.linspace(-dY/ddY,dY/ddY,nmax)
	xx, yy = np.meshgrid(x,y)

	maxradius = ( (dX/ddX)**2 + (dY/ddY)**2 )**0.5
	mult = 1                        # Increases or decreases the numbers of bins used. Most accurate results at mult=1.
	reselements = math.floor(mult*maxradius)
		                        # This is the number of "resolution elements" (i.e. the number of points
		                        #      on the corr_funct vs. radius plot) that we're dealing with.
	radiusmap = (xx**2 + yy**2)**0.5

	corr_functM51 = np.arange(nmax*reselements).reshape(nmax,reselements)
	corr_functM33 = np.arange(nmax*reselements).reshape(nmax,reselements)


	for i in range (0, dV/ddV+1):	# "i" defined as "dv/ddV".
		corr_functM51[i], edges, counts = ss.binned_statistic(
		radiusmap[radiusmap<maxradius], xi_M51[i][radiusmap<maxradius], statistic=np.nanmean, bins = reselements)
	for i in range (0, dV/ddV+1):	# "i" defined as "dv/ddV".
		corr_functM33[i], edges, counts = ss.binned_statistic(
		radiusmap[radiusmap<maxradius], xi_M33[i][radiusmap<maxradius], statistic=np.nanmean, bins = reselements)

	for j in range (0,np.int(reselements)):
		if np.isnan(corr_functM51[0,(reselements-1)-j]) == False:
			lowestM51 = corr_functM51[0,(reselements-1)-j]
		if np.isnan(corr_functM51[0,j]) == False:
			highestM51 = corr_functM51[0,j]
		if np.isnan(corr_functM33[0,(reselements-1)-j]) == False:
			lowestM33 = corr_functM33[0,(reselements-1)-j]
		if np.isnan(corr_functM33[0,j]) == False:
			highestM33 = corr_functM33[0,j]

	multM51 = highestM51
	multM33 = highestM33

	plt.figure(3)
	fig = matplotlib.pyplot.gcf()	
	fig.set_size_inches(15, 7)	# Enlarges the image so as to prevent squishing.

	X_M51 = (np.arange(reselements)/mult) / ((reselements-1)/mult) * (dX**2 + dY**2)**0.5 * pixelwidthPC_M51
	X_M33 = (np.arange(reselements)/mult) / ((reselements-1)/mult) * (dX**2 + dY**2)**0.5 * pixelwidthPC_M33
	for i in range (0, dV/ddV+1):							# Plot for M51.
	    if velocityresM51 > 0:
		plt.plot(X_M51, corr_functM51[i]/multM51,label='xi at +'+str('{0:.2f}'.format(i*ddV*velocityresM51))+' km/s for M51')
	    else:
		plt.plot(X_M51, corr_functM51[i]/multM51,label='xi at '+str('{0:.2f}'.format(i*ddV*velocityresM51))+' km/s for M51')

	for i in range (0, dV/ddV+1):							# Plot for M33.
	    if velocityresM33 > 0:
		plt.plot(X_M33, corr_functM33[i]/multM33,label='xi at +'+str('{0:.2f}'.format(i*ddV*velocityresM33))+' km/s for M33')
	    else:
		plt.plot(X_M33, corr_functM33[i]/multM33,label='xi at '+str('{0:.2f}'.format(i*ddV*velocityresM33))+' km/s for M33')

	print lowestM51
	print highestM51
	print lowestM33
	print highestM33

	plt.title('Avg. Corr. Funct. vs. Radial "Distance" from Center of xi Plots')
	plt.xlabel('Distance from Initial Location (pc)')
	plt.ylabel('Average xi (Normalized)')
	plt.legend(loc='upper left')
	plt.savefig('plot_M51andM33_comparison.png')
	plt.clf()


