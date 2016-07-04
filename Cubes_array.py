
# 6.08.16 - Works with S_2 surface maps that we've produced.


print('\nWelcome to Cubes_array! \n \nAvailable functions: \n  generate: Activates "anglefinder" and "slicer" for a given map.\n  slicer: Generates a thin "slice" of a S_2 surface map. \n  anglefinder: Finds the position angle of a given matrix. \n  plot: Generates a map showing the line of minimal S_2, along with a plot \n        showing the values along this line. \n  presetM51: Runs "generate" for various preset regions of M51.\n  presetM33: Same as above, but for M33.' )

from spectral_cube import SpectralCube
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tempfile import TemporaryFile
import scipy.interpolate as si
from scipy.signal import argrelextrema as ae

def generate(galaxyname='M51',vmin=40, vmax=80, ymin=200, ymax=400, xmin=360, xmax=560, deltaX=40, deltaV=3, deltadeltaX=1, deltadeltaV=1, nmax=201):
	"""
	Takes a S_2 surface map whose name matches the given parameters, finds the angle
		at which a line cutting through the origin has minimal S_2, and then
		returns a 1D array containing the (interpolated) values of S_2 along this
		line.

	Parameters:
	-----------
	galaxyname : string
		Name of the galaxy we're dealing with ('M51' or 'M33').
	vmin,...,deltadeltaV : int
		Parameters used in relevant S_2 map.
	nmax : int
		Number of values in final 1D array (preferably odd). 
		Higher is better.
	
	Returns (DISABLED at the moment):
	-----------
	theta : float
		Angle from x-axis at which a line cutting through 
		origin has minimal S_2.
	linearray1 : array (1D, float)	
		A list of S_2 values along this line, with the 
		origin in the middle of the list.
	linearrayx : array (1D, float)
		A list of "radius" values along this line, 
		corresponding to the data appearing in 'linearray1'.
		(It goes from -maxradius to +maxradius, where 
	  	'maxradius' is the maximum radius of the S_2 
		surface map.)
	"""
	imagename = galaxyname+"_"+str(vmin)+"to"+str(vmax)+"_"+str(ymin)+"to"+str(ymax)+"_"+str(xmin)+"to"+str(xmax)
	if deltadeltaX == 1 and deltadeltaV == 1:
		tempname = 'saved_S2array_'+imagename+'_dV_is_'+str(deltaV)+'_dX_is_'+str(deltaX)+'_MAXRES'
	else:
		tempname = 'saved_S2array_'+imagename+'_dV_is_'+str(deltaV)+'_dX_is_'+str(deltaX)

	f = file(tempname+".bin","rb")
	array = np.load(f)
	f.close()

	if galaxyname=='M51':
		filename = "paws_norot"
	else:
		if galaxyname=='M33':
			filename = "m33.co21_iram_CLEANED"
		else:
			print 'ERROR: Galaxy must be M51 or M33.'
			return

	theta = anglefinder(array[0].max()-array[0],False)
	linearrayx, linearray1, linearray2, maxradius1, maxradius2 = slicer(theta, array[0], nmax)
	plot(theta,maxradius1,maxradius2,array,linearrayx,linearray1,linearray2, filename, imagename, deltaX, deltaV)

#	return(theta,linearray1,linearrayx)

def anglefinder(weight, ReturnSizes=False):

	'''
	Calculates the position angle (measured counter-clockwise from the
	x-axis) for a weight matrix "weight".  This is the
	moment-of-inertia approach.
	
	Parameters:
	----------
	weight : float
		2D matrix showing weight for each point.  
		Should be non-negative.
	ReturnSizes : bool
		If True, return major and minor axis sizes along with 
		position angle.

	Returns:
	pa : float
		Position angle of the provided matrix.
	-------
	'''

	#
	#                _     _
	#               |       |
	#               | A   B |
	#      matrix = |       |
	#               | C   D |
	#               |_     _|
	#
	# where:
	#    A = sum of m_{ij}x_{ij}**2
	#    B = sum of m_{ij}y_{ij}x_{ij}
	#    C = B
	#    D = sum of m_{ij}y_{ij}**2
	# .    

	sumwts = np.nansum(weight)
	jmax,imax = weight.shape

	a = np.zeros(jmax*imax).reshape(jmax,imax)
	b = np.zeros(jmax*imax).reshape(jmax,imax)
	d = np.zeros(jmax*imax).reshape(jmax,imax)

	icen = (imax-1)/2.   # Central i-value, or central x-value.
	jcen = (jmax-1)/2.   # Central j-value, or central y-value.

	rmax = 0.35*np.min([jmax,imax])

	for j in range(0,jmax):
		for i in range(0,imax):
			if np.sqrt((i-icen)**2 + (j-jcen)**2) < rmax:
				a[j,i] = weight[j,i]*(i-icen)**2
				b[j,i] = weight[j,i]*(j-jcen)*(i-icen)
				d[j,i] = weight[j,i]*(j-jcen)**2
	A = np.nansum(a)
	B = np.nansum(b)
	C = B
	D = np.nansum(d)

	matrix = 1/sumwts*np.array([[A,B],[C,D]])
	determ = np.linalg.det(matrix)

	if ~np.isfinite(determ) or determ == 0:
		return(np.nan)
	evals, evecs = np.linalg.eigh(matrix)
	bigvec = evecs[-1,:]
	pa = np.arctan2(bigvec[1],bigvec[0])
	if ReturnSizes:
		major = evals[-1]
		minor = evals[0]
		return (pa,major,minor)

	return(pa)



def slicer(theta, array, nmax=201):
	"""
	Takes a S_2 surface map, cuts a thin "slice" through the middle at angle 
		'theta', and returns a 1D array containing the values of that slice.

	Parameters:
	-----------
	theta : float
		Angle from origin at which slicing should occur (rad).
	array : array
		The S_2 map that we're dealing with.
	nmax (int)
		Number of values in final 1D array. Higher
		is better.
	----------
	
	Returns:
	-----------
	linearrayx : array (1D, float)
		A list of "radius" values that corresponds
		to linearray1 and linearray2's data. It goes
		from -maxradius to +maxradius, where 'maxradius'
		is the maximum radius of the S_2 surface map.
	linearray1 : array (1D, float)	
		A list of S_2 values along the principal axis
		(i.e. the line of minimal S_2), with the 
		origin in the middle of the list.
	linearray2 : array (1D, float)
		Same as linearray1, but for the line
		perpendicular to the principal axis.
	maxradius1 : float
		The length of a line drawn along the principal
		axis from the center of the S_2 surface map
		to the map's outer edge. Important for plotting.
	maxradius2 : float
		Same as maxradius1, but for the line
		perpendicular to the principal axis.
	"""

	jmax, imax = array.shape

	X = np.arange(imax) - (imax-1)/2   # This is a vector going from (-deltaX,-deltaX+1,...,0,...,deltaX-1,deltaX).
	Y = np.arange(jmax) - (jmax-1)/2
	fxy = array                        # This is f(x,y). We'll interpolate so that it's a smooth surface,
		                           #        rather than being made up of countless "pixels".

	fxy1 = si.interp2d(X,Y,fxy)	   # This is f(x,y), or "array", but interpolated.
	#print fxy1(40,0)                  # Note the parentheses, not brackets.
		                           #  ALSO note: the coordinates are reversed!

	maxradius = np.sqrt( ((imax-1)/2)**2 + ((jmax-1)/2)**2 )    	# Largest "distance" from center of 'fxy'.
	maxradius1 = np.nan
	maxradius2 = np.nan

	nmax = nmax                                                 	# Must be odd for 'linearray1' to be perfectly 
		                                                    	#    zero in the middle. Other than that, even is fine.

	linearray1 = np.linspace(0,0,nmax)				# This will be the final array returned. It should be plotted against 'linearrayx', which is below.
	linearray2 = np.linspace(0,0,nmax)

	theta2 = theta + np.pi/2.

	for i in range(0,nmax):
	    r = (i-(nmax-1.)/2.) / ((nmax-1.)/2.) * maxradius
	    x = r*np.cos(theta)
	    y = r*np.sin(theta)
	    if np.abs(x) <= (imax-1.)/2. and np.abs(y) <= (jmax-1.)/2.:
		linearray1[i] = fxy1(x,y)
		maxradius1 = r                                      # Largest "distance" from center of 'fxy1' along the
		                                                    #    straight line at angle 'theta' from origin.
	    else:
		linearray1[i] = np.nan

	    x2 = r*np.cos(theta2)
	    y2 = r*np.sin(theta2)
	    if np.abs(x2) <= (imax-1.)/2. and np.abs(y2) <= (jmax-1.)/2.:
		linearray2[i] = fxy1(x2,y2)
		maxradius2 = r                                      # Largest "distance" from center of 'fxy1' along the
		                                                    #    straight line at angle 'theta' from origin PLUS 90 DEGREES.
	    else:
		linearray2[i] = np.nan
		
	linearrayx = np.linspace(-1,1,nmax) * maxradius


	return linearrayx,linearray1,linearray2,maxradius1,maxradius2


def plot(theta,maxradius1,maxradius2,array,linearrayx,linearray1,linearray2, filename="paws_norot", imagename="defaultname", deltaX=40, deltaV=3):
	'''
	Given a map, its position angle, the values along
		the low-S_2 "slice", and other information,
		generates a map showing the position angle
		and a plot showing the values along the 
		low-S_2 "slice".

	Parameters:
	-----------
	theta : float
		Angle of minimal S_2 along the provided map,
		counterclockwise from x-axis, in radians.
		Found using "Cubes_array.slicer".
	maxradius1 : float
		Distance between a line drawn at the map's
		centre (at angle 'theta') and the border of
		the map, in pixels.
		Found using "Cubes_array.slicer".
	maxradius2 : float
		Same as above, but at angle 'theta + pi/2'.
	array : array
		The S_2 map that we're dealing with.
	linearrayx : array (1D, float)
		A list of "radius" values along this line, 
		corresponding to the data appearing in 'linearray'.
		(It goes from -maxradius to +maxradius, where 
	  	'maxradius' is the maximum radius of the S_2 
		surface map.)
		Found using "Cubes_array.slicer".
	linearray1 : array (1D, float)
		A list of S_2 values along the line of
		minimal S_2, with the origin in the middle 
		of the list.
		Found using "Cubes_array.slicer".
	linearray2 : array (1D, float)
		A list of S_2 values perpendicular to the line of
		minimal S_2, with the origin in the middle 
		of the list.
		Found using "Cubes_array.slicer".
	filename : string
		Name of the .fits file that the array came
		from. "paws_norot" for M51, "m33.co21_iram_CLEANED"
		for M33.
	imagename : string
		Name of the image that will be saved.
	deltaX, deltaV : int
		Parameters used in the relevant S_2 map.
	'''

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.
	yshape = data.shape[1]/2.0
	xshape = data.shape[2]/2.0

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	if filename =='m33.co21_iram_CLEANED':			# Checks if the galaxy's Header file contains its distance.
	    distancePC = 840000.0				# The distance to the galaxy that M33's .fits file deals with, in parsecs. ONLY works on the CLEANED file!
	else:
	    distancePC = cube.header['DIST']		# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dV = deltaV		     	# Same as above, but for "dv".

	jmax,imax = array[0].shape

	theta2 = theta + np.pi/2.

	nmax = linearrayx.shape[0]

	dlinearray1 = np.gradient(linearray1, 2*maxradius1 / np.count_nonzero(~np.isnan(linearray1)) ,edge_order=2)		# Derivative of linearray1.
	ddlinearray1 = np.gradient(dlinearray1, 2*maxradius1 / np.count_nonzero(~np.isnan(linearray1)) ,edge_order=2)		# Second derivative of linearray1.

	linearray1_min = np.linspace(np.nan,np.nan,nmax)
	linearray1_min[ ae(ddlinearray1,np.greater,order=5) ] = linearray1[ ae(ddlinearray1,np.greater,order=5) ]				# Maxima of second derivative of linearray1.

	dlinearray2 = np.gradient(linearray2, 2*maxradius2 / np.count_nonzero(~np.isnan(linearray2)) ,edge_order=2)		# Same as above, but for linearray2.
	ddlinearray2 = np.gradient(dlinearray2, 2*maxradius2 / np.count_nonzero(~np.isnan(linearray2)) ,edge_order=2)		#

	linearray2_min = np.linspace(np.nan,np.nan,nmax)
	linearray2_min[ ae(ddlinearray2,np.greater,order=5) ] = linearray2[ ae(ddlinearray2,np.greater,order=5) ]				#



	# PLOTTING EVERYTHING
	# -------------------
	fig, axarr = plt.subplots(nrows=1,ncols=2)
	ax1, ax2 = axarr
	fig = plt.gcf()
	fig.set_size_inches(15,7.5)	# Enlarges the image so as to prevent squishing.

	ax1.imshow(array[0], interpolation = 'none', extent = [-dX*pixelwidthPC,dX*pixelwidthPC,-dY*pixelwidthPC,dY*pixelwidthPC],\
		   vmin=0, vmax=array.max(), aspect='auto', origin='lower')
	ax1.set_title('S_2 at 0 km/s')
	ax1.set_xlabel('Distance from Initial Location in x-direction (pc)')
	ax1.set_ylabel('Distance from Initial Location in y-direction (pc)')
	xmin = max(-maxradius1*np.cos(theta) , -(imax-1)/2 )
	xmax = min( maxradius1*np.cos(theta) , (imax-1)/2 )
	ymin = max(-maxradius1*np.sin(theta) , -(jmax-1)/2 )
	ymax = min( maxradius1*np.sin(theta) , (jmax-1)/2 )
	ax1.plot([xmin*pixelwidthPC,xmax*pixelwidthPC], [ymin*pixelwidthPC,ymax*pixelwidthPC], 'k-')
	xmin2 = max(-maxradius2*np.cos(theta2) , -(imax-1)/2 )
	xmax2 = min( maxradius2*np.cos(theta2) , (imax-1)/2 )
	ymin2 = max(-maxradius2*np.sin(theta2) , -(jmax-1)/2 )
	ymax2 = min( maxradius2*np.sin(theta2) , (jmax-1)/2 )
	ax1.plot([xmin2*pixelwidthPC,xmax2*pixelwidthPC], [ymin2*pixelwidthPC,ymax2*pixelwidthPC], 'k:')

	ax1.set_xlim(-dX*pixelwidthPC-0.5,dX*pixelwidthPC-0.5)
	ax1.set_ylim(-dY*pixelwidthPC+0.5,dY*pixelwidthPC+0.5)


	ax2.plot(linearrayx*pixelwidthPC,linearray1,'k-',label='Principal Axis')
	ax2.plot(linearrayx*pixelwidthPC,linearray1_min,'r.')
	ax2.plot(linearrayx*pixelwidthPC,linearray2,'k:',label='Principal Axis + pi/2')
	ax2.plot(linearrayx*pixelwidthPC,linearray2_min,'g.')
	ax2.set_title('S_2 along "Line of Lowest S_2", versus Radius')
	ax2.set_xlabel('Distance from Center along Line (pc)')
	ax2.set_ylabel('S_2')
	ax2.legend(loc='lower left')

	plt.tight_layout()
	plt.savefig("S2_minimal_"+imagename+".png")
	plt.clf()
	# -------------------

def presetM51(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=1, deltadeltaV=1):
	'''
	Activates Cubes_array.generate using each of
		the preset regions of M51.

	Parameters:
	-----------
	vmin,...,deltadeltaV : int
		Parameters used in the relevant S_2 map.
	'''
	
	galaxyname = 'M51'

	ymin = np.array([350,200,220,350,350,100,200])	# These are the minimum "y" values of the regions that we're dealing with.
	ymax = np.array([550,400,420,550,550,300,400])	# These are the corresponding maximum "y" values of these regions.
	xmin = np.array([500,425,260,120,250,570,360])	# These are the corresponding minimum "x" values of these regions.
	xmax = np.array([700,625,460,320,450,770,560])	# These are the corresponding maximum "x" values of these regions. (Example: The first region has ymin=350, ymax=550, xmin=500, xmax=700.)
	sets = np.ravel(ymin.shape)[0]		# This is the number of regions that we're dealing with.

	for i in range(0,sets):
		generate(galaxyname, vmin,vmax,ymin[i],ymax[i],xmin[i],xmax[i],deltaX,deltaV,deltadeltaX,deltadeltaV,201)


def presetM33(vmin=40,vmax=80, deltaX=40, deltaV=6, deltadeltaX=1, deltadeltaV=1):
	'''
	Activates Cubes_array.generate using each of
		the preset regions of M33.

	Parameters:
	-----------
	vmin,...,deltadeltaV : int
		Parameters used in the relevant S_2 map.
	'''
	
	galaxyname = 'M33'

	ymin = np.array([350,600,650,525,300,250])	# These are the minimum "y" values of the regions that we're dealing with.
	ymax = np.array([550,800,850,725,500,450])	# These are the corresponding maximum "y" values of these regions.
	xmin = np.array([500,100,400,288,200,550])	# These are the corresponding minimum "x" values of these regions.
	xmax = np.array([700,300,600,488,400,750])	# These are the corresponding maximum "x" values of these regions. (Example: The first region has ymin=350, ymax=550, xmin=500, xmax=700.)
	sets = np.ravel(ymin.shape)[0]		# This is the number of regions that we're dealing with.

	for i in range(0,sets):
		generate(galaxyname, vmin,vmax,ymin[i],ymax[i],xmin[i],xmax[i],deltaX,deltaV,deltadeltaX,deltadeltaV,201)

