
# 6.08.16 - Works with S_2 surface maps that we've produced.


print('\nWelcome to Cubes_array! \n \nAvailable functions: \n  slicer: Generates a thin "slice" of a S_2 surface map. \n  anglefinder: Finds the position angle of a given matrix.' )

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from tempfile import TemporaryFile
import scipy.interpolate as si

def anglefinder(weight, ReturnSizes=False):
	sumwts = np.nansum(weight)

	'''
	Calculates the position angle (measured counter-clockwise from the
	x-axis) for a weight matrix "weight".  This is the
	moment-of-inertia approach.
	
	Parameters
	----------
	weight : float
		2D matrix showing weight for each point.  
		Should be non-negative.
	ReturnSizes : bool
		If True, return major and minor axis sizes along with 
		position angle.
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

	jmax,imax = weight.shape

	a = np.zeros(jmax*imax).reshape(jmax,imax)
	b = np.zeros(jmax*imax).reshape(jmax,imax)
	d = np.zeros(jmax*imax).reshape(jmax,imax)

	icen = (imax-1)/2.   # Central i-value, or central x-value.
	jcen = (jmax-1)/2.   # Central j-value, or central y-value.

	for j in range(0,jmax):
	for i in range(0,imax):
	    a[j,i] = weight[j,i]*(i-icen)**2
	    b[j,i] = weight[j,i]*(j-jcen)*(i-icen)
	    d[j,i] = weight[j,i]*(j-jcen)**2

	A = np.nanmean(a)
	B = np.nanmean(b)
	C = B
	D = np.nanmean(d)

	matrix = 1/sumwts*np.array([[A,B],
		                [C,D]])

	plt.imshow(matrix,interpolation='none')
	plt.colorbar()

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
	return(np.pi-pa)



def slicer(theta, galaxyname='M51',vmin=40, vmax=80, ymin=200, ymax=400, xmin=360, xmax=560, deltaX=40, deltaV=3, deltadeltaX=1, deltadeltaV=1, nmax=201):
	"""
	Takes a S_2 surface map whose name matches the given parameters, cuts a thin "slice" through the middle at
		angle 'theta', and returns a 1D array containing the values of that slice.

	Parameters:
		theta 		- Angle from origin at which slicing should occur (rad).
		galaxyname	- Name of the galaxy we're dealing with.
		vmin,vmax,ymin,ymax,xmin,xmax,deltaX,deltaV,deltadeltaX,deltadeltaV
				- Parameters used in relevant S_2 map.
		nmax		- Number of values in final 1D array. Higher is better.
	"""

	imagename = galaxyname+"_"+str(vmin)+"to"+str(vmax)+"_"+str(ymin)+"to"+str(ymax)+"_"+str(xmin)+"to"+str(xmax)
	if deltadeltaX == 1 and deltadeltaV == 1:
		tempname = 'saved_S2array_'+imagename+'_dV_is_'+str(deltaV)+'_dX_is_'+str(deltaX)+'_MAXRES'
	else:
		tempname = 'saved_S2array_'+imagename+'_dV_is_'+str(deltaV)+'_dX_is_'+str(deltaX)

	

	f = file(tempname+".bin","rb")
	array = np.load(f)
	f.close()

	jmax, imax = array[0].shape

	X = np.arange(imax) - (imax-1)/2   # This is a vector going from (-deltaX,-deltaX+1,...,0,...,deltaX-1,deltaX).
	Y = np.arange(jmax) - (jmax-1)/2
	fxy = array                        # This is f(x,y). We'll interpolate so that it's a smooth surface,
		                           #        rather than being made up of countless "pixels".

	fxy1 = si.interp2d(X,Y,array[0])   # This is f(x,y), or "array", but interpolated.
	fxy1(40,0)                         # Note the parentheses, not brackets.
		                           #  ALSO note: the coordinates are reversed!

	maxradius = np.sqrt( ((imax-1)/2)**2 + ((jmax-1)/2)**2 )    	# Largest "distance" from center of 'fxy'.
	nmax = nmax                                                 	# Must be odd for 'linearray' to be perfectly 
		                                                    	#    zero in the middle. Other than that, even is fine.

	linearray = np.linspace(0,0,nmax)				# This will be the final array returned. It should be plotted against 'linearrayx', which is below.

	for i in range(0,nmax):
	    r = (i-(nmax-1.)/2.) / ((nmax-1.)/2.) * maxradius
	    x = r*np.cos(theta)
	    y = r*np.sin(theta)
	    if np.abs(x) <= (imax-1.)/2. and np.abs(y) <= (jmax-1.)/2.:
		linearray[i] = fxy1(x,y)
		maxradius2 = r                                      # Largest "distance" from center of 'fxy1' along the
		                                                    #    straight line at angle 'theta' from origin.
	    else:
		linearray[i] = np.nan
		
	linearrayx = np.linspace(-1,1,nmax) * maxradius
#	plt.plot(linearrayx,linearray,'.')
#	plt.show()
#	plt.clf()
	return linearray#,linearrayx
