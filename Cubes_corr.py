
# 6.15.16 - Same as "Cubes", but deals with the CORRELATION function instead of the Structure function.

from spectral_cube import SpectralCube
import matplotlib
import matplotlib.pyplot as plt
import astropy.units as u
import numpy as np
import scipy.stats as ss
import math

print('"Cubes_corr.py" \n \nAvailable functions within Cubes_corr: \n  cubegen - generates a subcube. \n  corrgen - generates a correlation function map from a given subcube. \n  mapgen - generates some 2D maps of the correlation function, for each "dv". \n  plotgen - generates a 1D plot of averaged correlation function versus radius. \n  everythinggen - generates 3-panel image with Tmax map, xi map, 1D plot.')

def cubegen(vmin,vmax,ymin,ymax,xmin,xmax, filename = "paws_norot", drawmap = False, mapname="3Dcube"):
	"""
	Returns a subcube of the specified dimensions from the .fits file.
	Also displays the subcube as it appears on the galaxy map if drawmap=True.


	Parameters:
	-----------
	vmin,...,xmax : int
		Parameters used in relevant xi map.
		WARNING: Selecting too large of a vmax-vmin will hugely increase
		processing time in later calculations.
	filename : str
		Name of the .paws data file.
		"paws_norot" for M51, "m33.co21_iram_CLEANED" for M33.
	drawmap : bool
		Enables or disables drawing the subcube Tmax map.
	galaxyname : str
		Name of the galaxy.
		'M51' for M51, 'M33' for M33.
	mapname : str
		Name of the saved image of the subcube's Tmax map, if
		drawmap==True.


	Returns:
	-----------
	subcube : spectral cube (?)
		The data inside the selected subcube.
"""

	cube = SpectralCube.read(filename+".fits")
	data = cube.filled_data[:]   # Pulls "cube"'s information (position, spectral info (?)) into a 3D Numpy array.
	yshape = data.shape[1]/2.0
	xshape = data.shape[2]/2.0

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	if (filename =='m33.co21_iram_CLEANED') or (filename =='m33.co21_iram_CLEANED_smooth') or (filename =='m33.co21_iram_CLEANED_blank'):	# Checks if the galaxy's Header file contains its distance.
		distancePC = 840000.0				# The distance to the galaxy that M33's .fits file deals with, in parsecs. ONLY works on the CLEANED file!
	else:
		distancePC = cube.header['DIST']		# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	subcube = cube[vmin:vmax,ymin:ymax,xmin:xmax]
	if drawmap == True:
		plt.figure(1)
		plt.imshow(np.nanmax(data[vmin:vmax,ymin:ymax,xmin:xmax].value,axis=0), extent=[(xmin-xshape)*pixelwidthPC,(xmax-xshape)*pixelwidthPC, \
			   (ymin-yshape)*pixelwidthPC,(ymax-yshape)*pixelwidthPC], origin='lower')
		fig = matplotlib.pyplot.gcf()
		#fig.set_size_inches(5, 5)	# Enlarges the image so as to prevent squishing.
		plt.xlabel('Distance from Centre in x-direction (pc)')
		plt.ylabel('Distance from Centre in y-direction (pc)')

		plt.savefig('galaxy_'+mapname+'.png')
		plt.clf()			# Clears the image after saving.

	return subcube

def corrgen(subcube0, deltaX=30, deltaV=3, deltadeltaX=1, deltadeltaV = 1):
	"""
	Returns a 3D correlation function map from a given subcube.


	Parameters:
	-----------
	subcube : spectral cube (?)
		The subcube generated in cubegen.
	deltaX,...,deltadeltaV : int
		Parameters used in relevant xi map.
		WARNING: Selecting too large of a vmax-vmin will hugely increase
		processing time.

	Returns:
	-----------
	xi : array
		3D map of the correlation function. The third dimension
		represents shift in the spectral axis."""

	data0 = subcube0.filled_data[:]

	dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dV = deltaV		     	# Same as above, but for "dv".
	ddX = deltadeltaX		# "Step size" (in units of spatial resolutions), for the xi calculation. Make sure it's a factor of dX. Leave at 1 for full calculation.
	ddY = np.copy(ddX)            	# Same as above, but for dY.
	ddV = deltadeltaV		# Same as above, but for dV (in units of spectral resolutions). Not the same as ddX or ddY.
	nmax = abs(2*dX/ddX)+1
	xi = np.zeros([np.abs(dV/ddV)+1,nmax,nmax])

	p,n,m = data0.shape     # The subcube has n rows, m columns, and p "height units".


	if dV>=0:
		for delv in range (0,dV/ddV+1):                 # "delv" defined as "dv / ddV"
		    for delx in range (-dX/ddX,dX/ddX+1):       # "delx" defined as "dx / ddX"
			for dely in range (-dY/ddY,dY/ddY+1):   # "dely" defined as "dy / ddY"

			    dx = delx*ddX
			    dy = dely*ddY
			    dv = delv*ddV
		
			    M = data0.value         		# This is the subcube information (unitless).
			    P = np.arange(p*n*m).reshape(p,n,m) # This will be used to track the shifting "pixels" of M(r) and M(r+dr).
			    D = np.zeros([p,n,m])   		# This will be the difference between M(r) and M(r+dr).

			    RollMap = np.roll((np.roll(np.roll(M,-dv,axis=0),-dy,axis=1)),-dx,axis=2)
			    D = M - RollMap

			    goodpix = (P - np.roll(P,-dv,axis=0) == -dv*n*m) * (P - np.roll(P,-dy,axis=1) == -dy*m) * (P - np.roll(P,-dx,axis=2) == -dx)
				    # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
				    #        pixel above or below it by exactly m. So, the difference between a pixel's value and
				    #        that of a pixel "dy" rows below is simply dy*m.
				    #	    Similar case for a pixel "dv" layers below.
				    #
				    # In "goodpix", pixels that have wrapped around are treated as False.

	#		    OrigMapPower = np.nanmean(M[goodpix]**2)
	#		    RollMapPower = np.nanmean(RollMap[goodpix]**2)

			    xi[(dv/ddV,(dy+dY)/ddY,(dx+dX)/ddX)] = (np.nanmean(D[goodpix]**2) - np.nanmean(M**2) - np.nanmean(RollMap**2)) / (-2*np.nanmean(M)**2)
	#		    Power = (OrigMapPower + RollMapPower)
	else:	# If dV < 0:
		for delv in range (dV/ddV,0+1,1):                 	# "delv" defined as "dv / ddV"		# <-- (!) SAME AS ABOVE, but in the NEGATIVE "direction".
		    for delx in range (-dX/ddX,dX/ddX+1):       	# "delx" defined as "dx / ddX"
			for dely in range (-dY/ddY,dY/ddY+1):   	# "dely" defined as "dy / ddY"

			    dx = delx*ddX
			    dy = dely*ddY
			    dv = delv*ddV
		
			    M = data0.value         		# This is the subcube information (unitless).
			    P = np.arange(p*n*m).reshape(p,n,m) # This will be used to track the shifting "pixels" of M(r) and M(r+dr).
			    D = np.zeros([p,n,m])   		# This will be the difference between M(r) and M(r+dr).

			    RollMap = np.roll((np.roll(np.roll(M,-dv,axis=0),-dy,axis=1)),-dx,axis=2)
			    D = M - RollMap

			    goodpix = (P - np.roll(P,-dv,axis=0) == -dv*n*m) * (P - np.roll(P,-dy,axis=1) == -dy*m) * (P - np.roll(P,-dx,axis=2) == -dx)
				    # Note: The "-dy*m" is because, for P, a pixel's value is separated from that of a
				    #        pixel above or below it by exactly m. So, the difference between a pixel's value and
				    #        that of a pixel "dy" rows below is simply dy*m.
				    #	    Similar case for a pixel "dv" layers below.
				    #
				    # In "goodpix", pixels that have wrapped around are treated as False.

	#		    OrigMapPower = np.nanmean(M[goodpix]**2)
	#		    RollMapPower = np.nanmean(RollMap[goodpix]**2)

			    xi[(-dv/ddV,(dy+dY)/ddY,(dx+dX)/ddX)] = (np.nanmean(D[goodpix]**2) - np.nanmean(M**2) - np.nanmean(RollMap**2)) / (-2*np.nanmean(M)**2)	# <-- (!) "-dv/ddV" so that xi[0]
																					#    is dv=0, and xi[|dV|] is dv=dV.
	#		    Power = (OrigMapPower + RollMapPower)

	return xi


def mapgen(xi, deltaX=30, deltaV=3, deltadeltaV=1, mapname="3Dcube", filename = "paws_norot"):
	"""
	Generates and saves several 2D colour plots of the correlation function versus
	   position; one for no shift in spectral dimension and one for "maximum" shift
	   in spectral dimension.

	Parameters:
	-----------
	xi : array
		3D map of the correlation function, generated in corrgen.
	deltaX,...,deltadeltaV : int
		Parameters used in relevant xi map.
	mapname : str
		Name of the saved image of the subcube's xi map.
	filename : str
		Name of the .paws data file.
		"paws_norot" for M51, "m33.co21_iram_CLEANED" for M33.

	Returns:
	-----------
	none
	"""

	dX = deltaX                  	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dV = deltaV		     	# Same as above, but for "dv".
	ddV = deltadeltaV
	cube = SpectralCube.read(filename+".fits")

	pixelwidthDEG = cube.header['CDELT2']			# The width of each pixel, in degrees.
	if filename =='m33.co21_iram_CLEANED':			# Checks if the galaxy's Header file contains its distance.
		distancePC = 840000.0				# The distance to the galaxy that M33's .fits file deals with, in parsecs. ONLY works on the CLEANED file!
	else:
		distancePC = cube.header['DIST']		# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
	pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.

	velocityres = cube.header['CDELT3']	# Velocity resolution in m/s or km/s.
	if filename != 'paws_norot':                 # 'paws_norot' already has its velocity resolution in km/s.
	    velocityres = velocityres / 1000.0
	# (!)					# May need to edit these above three lines for other .fits files!

	# 2D Display (layer by layer)
	plt.figure(2)

	fig = matplotlib.pyplot.gcf()
	fig.set_size_inches(15, 7)

	plt.subplot(121)
	plt.imshow(xi[0], interpolation = 'none', extent = [-dX*pixelwidthPC,dX*pixelwidthPC,-dY*pixelwidthPC,dY*pixelwidthPC], vmin=0, vmax=xi.max(), aspect='auto', origin='lower')
	plt.title('xi at 0 km/s')
	#plt.xlabel('Distance from Initial Location in x-direction (pc)')
	#plt.ylabel('Distance from Initial Location in y-direction (pc)')

	plt.subplot(122)
	plt.imshow(xi[np.abs(dV/ddV)], interpolation = 'none', extent = [-dX*pixelwidthPC,dX*pixelwidthPC,-dY*pixelwidthPC,dY*pixelwidthPC], vmin=0, vmax=xi.max(), aspect='auto', origin='lower')
	if dV*velocityres > 0:
	    plt.title('xi at +'+str('{0:.2f}'.format(dV*velocityres))+' km/s')
	else:
	    plt.title('xi at '+str('{0:.2f}'.format(dV*velocityres))+' km/s')
	#plt.xlabel('Distance from Initial Location in x-direction (pc)')
	#plt.ylabel('Distance from Initial Location in y-direction (pc)')

	plt.text(-dX*pixelwidthPC*1.2, -dY*pixelwidthPC*1.2, 'Distance from Initial Location in x-direction (pc)', ha='center', va='center')
	plt.text(-dX*pixelwidthPC*4.3, 0, 'Distance from Initial Location in y-direction (pc)', ha='center', va='center', rotation='vertical')

	plt.colorbar()
	plt.savefig('map_xi_'+mapname+'.png')
	plt.clf()			# Clears the figure, allowing "Figure 2" to be used again if a function calls on mapgen more than once.

def plotgen(xi, deltaX=30, deltaV=3, deltadeltaX=1, deltadeltaV=3, mapname="3Dcube", filename="paws_norot", J = 10000):
	"""
	Generates and saves a 1D plot of the azimuthally-averaged correlation function  
	   versus radius, for each value of "dv".

	Parameters:
	-----------
	xi : array
		3D map of the correlation function, generated in corrgen.
	deltaX,...,deltadeltaV : int
		Parameters used in relevant xi map.
		Setting deltaX=0 will use velocity shift instead of
		position shift.
	mapname : str
		Name of the saved image of the subcube's xi plot.
	filename : str
		Name of the .paws data file.
		"paws_norot" for M51, "m33.co21_iram_CLEANED" for M33.
	J : int
		Number of loops used in the bootstrap algorithm (which is
		used to generate error values for slope of log(xi) vs.
		log(shift).

	Returns:
	-----------
	coeff_a : float
		Intercept of log(xi) vs. log(shift).
	coeff_b : float
		Slope of log(xi) vs. log(shift).
	a_error : float
		Standard error of coeff_a.
	b_error : float
		Standard error of coeff_b.

	Note that this "shift" may be the shift of position OR radial
		velocity, depending on your setting for deltaX."""



	if deltaX != 0:
		# Goal: Create a 1D plot, for each value of "dv", of the average value of correlation function (inside a thin ring
		#       at radius r) versus radius. Each "dv" is a different sheet of dx,dy.

		dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
		dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
		dV = deltaV		     	# Same as above, but for "dv".
		ddX = deltadeltaX
		ddY = np.copy(ddX)
		ddV = deltadeltaV
		nmax = abs(2*dX/ddX)+1

		cube = SpectralCube.read(filename+".fits")
		pixelwidthDEG = cube.header['CDELT2']	# The width of each pixel, in degrees.
		if filename =='m33.co21_iram_CLEANED':			# Checks if the galaxy's Header file contains its distance.
			distancePC = 840000.0				# The distance to the galaxy that M33's .fits file deals with, in parsecs. ONLY works on the CLEANED file!
		else:
			distancePC = cube.header['DIST']		# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
		pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.
		velocityres = cube.header['CDELT3']	# Velocity resolution in m/s or km/s.
		if filename != 'paws_norot':            # 'paws_norot' already has its velocity resolution in km/s.
		    velocityres = velocityres / 1000.0
		# (!)					# May need to edit these above three lines for other .fits files!

		x = np.linspace(-dX/ddX,dX/ddX,nmax)
		y = np.linspace(-dY/ddY,dY/ddY,nmax)
		xx, yy = np.meshgrid(x,y)

		maxradius = ( (dX/ddX)**2 + (dY/ddY)**2 )**0.5
		mult = 1                        # Increases or decreases the numbers of bins used. Most accurate results at mult=1.
		reselements = math.floor(mult*maxradius)
				                # This is the number of "resolution elements" (i.e. the number of points
				                #      on the corr_funct vs. radius plot) that we're dealing with.
		radiusmap = (xx**2 + yy**2)**0.5

		corr_funct = np.arange(nmax*reselements).reshape(nmax,reselements)

		for i in range (0, np.abs(dV/ddV)+1):	# "delv" defined as "dv/ddV".
			corr_funct[i], edges, counts = ss.binned_statistic(
			radiusmap[radiusmap<maxradius], xi[i][radiusmap<maxradius], statistic=np.nanmean, bins = reselements)

		plt.figure(3)
		fig = matplotlib.pyplot.gcf()	
		fig.set_size_inches(15, 7)	# Enlarges the image so as to prevent squishing.
		x_axis = (np.arange(reselements)/mult) / ((reselements-1)/mult) * (dX**2 + dY**2)**0.5 * pixelwidthPC

		# Creates the plot.
		for i in range (0, np.abs(dV/ddV)+1):
		    if dV*velocityres > 0:
			plt.plot(x_axis[x_axis<dX*pixelwidthPC], corr_funct[i][x_axis<dX*pixelwidthPC],label='xi at +'+str('{0:.2f}'.format(i*ddV*velocityres*dV/np.abs(dV)))+' km/s')
		    elif dV*velocityres < 0:
			plt.plot(x_axis[x_axis<dX*pixelwidthPC], corr_funct[i][x_axis<dX*pixelwidthPC],label='xi at '+str('{0:.2f}'.format(i*ddV*velocityres*dV/np.abs(dV)))+' km/s')
		    else:
			plt.plot(x_axis[x_axis<dX*pixelwidthPC], corr_funct[i][x_axis<dX*pixelwidthPC],label='xi at '+str('{0:.2f}'.format(i*ddV*velocityres))+' km/s')
		plt.title('Avg. Corr. Funct. vs. Radial "Distance" from Center of xi Plots')
		plt.xlabel('Distance from Initial Location (pc)')
		plt.ylabel('Average xi')
	#	plt.ylim([0.9,1.1])

		# Turns the plot into log-log form as long as there are some positive values.
		positive_value = False				# Checks if there are positive xi values.
		for i in range (0, np.abs(dV/ddV)+1):
			if np.nanmax(corr_funct[i][x_axis<dX*pixelwidthPC]) > 0:
				positive_value = True

		if positive_value == False:
			print "~~~   ERROR: NO POSITIVE VALUES IN CORRELATION FUNCTION WHATSOEVER.  ~~~"
			plt.savefig('plot_xi_'+mapname+'.png')
			plt.clf()
			return np.nan, np.nan, np.nan, np.nan		
		else:
			plt.yscale('log')
			plt.xscale('log')

		# Calculates the intercept (a) and slope (b) of log(average "xi") vs. log(radial distance) at distances between 50pc and 250pc.
		Y = corr_funct[0][x_axis<dX*pixelwidthPC]
		X = x_axis[x_axis<dX*pixelwidthPC]

		Y = Y[X>=50]
		X = X[X>=50]
		Y = Y[X<=250]
		X = X[X<=250]
		X = np.log10(X)
		Y = np.log10(Y)
		gooddata = np.isfinite(X)*np.isfinite(Y)
		if (X[gooddata].size > 0) and (Y[gooddata].size > 0):
			coeff_b, coeff_a = np.polyfit(X[gooddata], Y[gooddata], 1)
		else:
			print "ERROR: X and/or Y do not have any 'good' values."
			return np.nan, np.nan, np.nan, np.nan
		print coeff_a, coeff_b

		# Generates standard deviation of slope (and, from there, error bars) using a bootstrap algorithm.
		b_trial = [None]*J		# "J" is the number of iterations the bootstrap algorithm will loop through.
		a_trial = [None]*J

		X_good = X[gooddata]
		Y_good = Y[gooddata]
		index = np.arange(X_good.size)

		for j in range(0,J):
		    choices = np.random.choice(index, size=index.size, replace=True)
		    b_trial[j],a_trial[j] = np.polyfit(X_good[choices],Y_good[choices],1)

		b_error = np.std(b_trial)
		a_error = np.std(a_trial)


		# Plots the line (or, without the log-log plot, curve) of best fit.
		if velocityres>0:
			plt.plot(10**X,10**(coeff_a + coeff_b*X),'k--',label='Best fit, at +'+str('{0:.2f}'.format(0*ddV*velocityres))+' km/s')
		else:
			plt.plot(10**X,10**(coeff_a + coeff_b*X),'k--',label='Best fit, at '+str('{0:.2f}'.format(0*ddV*velocityres))+' km/s')
		plt.legend(loc='upper right')

		if coeff_b > 0:
			plt.figtext(.15,.15,"Best fit, from 50pc to 250pc: \nlog(xi) = "+str(coeff_a)+" + "+str(coeff_b)+"log(radius)")
		else:
			plt.figtext(.15,.15,"Best fit, from 50pc to 250pc: \nlog(xi) = "+str(coeff_a)+" - "+str(np.abs(coeff_b))+"log(radius)")
		plt.savefig('plot_xi_'+mapname+'.png')
		plt.clf()

		# Plots and saves a histogram of the trial slopes from earlier.
#		print b_trial
		plt.hist(b_trial,bins=75)
		plt.title('Distribution of Trial Slopes')
		plt.xlabel('Trial Slope')
		plt.ylabel('# of Appearances')
		plt.savefig('plot_hist_xi_'+mapname+'.png')
		plt.clf()

		return coeff_a, coeff_b, a_error, b_error

	elif deltaX==0:
		# Goal: Create a 1D plot of the correlation function at each shift in velocity resolution element "dv". Basically, "xi" vs. "shift in radial velocity".

		dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
		dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
		dV = deltaV		     	# Same as above, but for "dv".
		ddV = deltadeltaV
		nmax = 1

		cube = SpectralCube.read(filename+".fits")
		pixelwidthDEG = cube.header['CDELT2']	# The width of each pixel, in degrees.
		if filename =='m33.co21_iram_CLEANED':			# Checks if the galaxy's Header file contains its distance.
			distancePC = 840000.0				# The distance to the galaxy that M33's .fits file deals with, in parsecs. ONLY works on the CLEANED file!
		else:
			distancePC = cube.header['DIST']		# The distance to the galaxy that M51's .fits file deals with, in parsecs.  (???) Is this number accurate, though?
		pixelwidthPC = pixelwidthDEG*np.pi/180.0*distancePC	# The width of each pixel, in pc.
		velocityres = cube.header['CDELT3']	# Velocity resolution in m/s or km/s.
		if filename != 'paws_norot':            # 'paws_norot' already has its velocity resolution in km/s.
		    velocityres = velocityres / 1000.0
		# (!)					# May need to edit these above three lines for other .fits files!

		corr_funct = xi
		corr_funct.shape = corr_funct.shape[0]

		print xi.shape

		plt.figure(3)
		fig = matplotlib.pyplot.gcf()	
		fig.set_size_inches(15, 7)			# Enlarges the image so as to prevent squishing.

		V = (np.arange(np.abs(dV/ddV)+1)) * velocityres	* (dV/np.abs(dV))	# This is an array of the different velocity shift values.
		if (V[0]<=0) and (V[-1]<0):
		    V = V*-1					# Flips sign of x-axis, so that the log of xi can be taken properly.

		    plt.plot(V, corr_funct,label='xi vs. NEGATIVE radial velocity shift')
		    plt.title('Corr. Funct. vs. "Radial Velocity Shift" from Center of xi Plots')
		    plt.xlabel('Shift in Radial Velocity in the NEGATIVE Direction (km/s)')
		    plt.ylabel('xi')
		else:
		    plt.plot(V, corr_funct,label='xi vs. radial velocity shift')
		    plt.title('Corr. Funct. vs. "Radial Velocity Shift" from Center of xi Plots')
		    plt.xlabel('Shift in Radial Velocity (km/s)')
		    plt.ylabel('xi')
		if np.nanmax(corr_funct)>0:
			plt.yscale('log')
			plt.xscale('log')
		else:
			print "ERROR!: "+mapname

		# Calculates the intercept (a) and slope (b) of log(average "xi") vs. log(Radial Velocity Shift).
		Y = corr_funct
		X = V

		X = np.log10(X)
		Y = np.log10(Y)


		gooddata = np.isfinite(X)*np.isfinite(Y)
		coeff_b, coeff_a = np.polyfit(X[gooddata], Y[gooddata], 1)

		# Generates standard deviation of slope (and, from there, error bars) using a bootstrap algorithm.
		b_trial = [None]*J		# "J" is the number of iterations the bootstrap algorithm will loop through.
		a_trial = [None]*J

		X_good = X[gooddata]
		Y_good = Y[gooddata]
		index = np.arange(X_good.size)

		for j in range(0,J):
		    choices = np.random.choice(index, size=index.size, replace=True)
		    b_trial[j],a_trial[j] = np.polyfit(X_good[choices],Y_good[choices],1)

		b_error = np.std(b_trial)
		a_error = np.std(a_trial)

		# Plots the line (or, without the log-log plot, curve) of best fit.

		plt.plot(10**X,10**(coeff_a + coeff_b*X),'k--',label='Best fit')
		plt.legend(loc='upper right')

		if coeff_b > 0:
			plt.figtext(.15,.15,"Best fit: \nlog(xi) = "+str(coeff_a)+" + "+str(coeff_b)+"log(shift)")
		else:
			plt.figtext(.15,.15,"Best fit: \nlog(xi) = "+str(coeff_a)+" - "+str(np.abs(coeff_b))+"log(shift)")
		plt.savefig('plot_xi_'+mapname+'.png')
		plt.clf()

		return coeff_a, coeff_b, a_error, b_error


def everythinggen(vmin, vmax, ymin, ymax, xmin, xmax, xi, deltaX, deltaV, deltadeltaX, deltadeltaV, imagename, filename):
	"""
	Generates and saves a 3-panel image showing a Tmax map of the
		selected region, a map of correlation function versus
		position shift, and a plot of correlation function versus
		position shift.

	Parameters:
	-----------
	vmin,...,deltadeltaV : int
		Parameters used in relevant xi map.
	imagename : str
		Name of the final saved image.
	filename : str
		Name of the .paws data file.
		"paws_norot" for M51, "m33.co21_iram_CLEANED" for M33.

	Returns:
	-----------
	none
	"""
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

	velocityres = cube.header['CDELT3']			# Velocity resolution in m/s or km/s.
	if filename != 'paws_norot':                 		# 'paws_norot' already has its velocity resolution in km/s.
	    velocityres = velocityres / 1000.0

	dX = deltaX                    	# This is simply the maximum absolute value of "dx". So if dX = 1, then dx = {-1,0,1}.
	dY = np.copy(dX)                # Same as above, but for "dy". For simplicity, let it be the same as dX.
	dV = deltaV		     	# Same as above, but for "dv".
	ddX = deltadeltaX
	ddY = np.copy(ddX)
	ddV = deltadeltaV
	nmax = abs(2*dX/ddX)+1

	x = np.linspace(-dX/ddX,dX/ddX,nmax)
	y = np.linspace(-dY/ddY,dY/ddY,nmax)
	xx, yy = np.meshgrid(x,y)

	maxradius = ( (dX/ddX)**2 + (dY/ddY)**2 )**0.5
	mult = 1                        # Increases or decreases the numbers of bins used. Most accurate results at mult=1.
	reselements = math.floor(mult*maxradius)
		                        # This is the number of "resolution elements" (i.e. the number of points
		                        #      on the corr_funct vs. radius plot) that we're dealing with.
	radiusmap = (xx**2 + yy**2)**0.5

	corr_funct = np.arange(nmax*reselements).reshape(nmax,reselements)

	for i in range (0, np.abs(dV/ddV)+1):	# "delv" defined as "dv/ddV".
		corr_funct[i], edges, counts = ss.binned_statistic(
		radiusmap[radiusmap<maxradius], xi[i][radiusmap<maxradius], statistic=np.nanmean, bins = reselements)
	X = (np.arange(reselements)/mult) / ((reselements-1)/mult) * (dX**2 + dY**2)**0.5 * pixelwidthPC


	fig, axarr = plt.subplots(nrows=1,ncols=3)
	#fig.subplots_adjust(wspace=0.8,bottom=0.2,right=0.85)
	#fig.set_figsize=(10,4)
	ax1, ax2, ax3 = axarr
	fig = plt.gcf()
	fig.set_size_inches(20,6)	# Enlarges the image so as to prevent squishing.

	#ax1 = fig.add_subplot(131)
	### Map
	ax1.imshow(np.nanmax(data[vmin:vmax,ymin:ymax,xmin:xmax].value,axis=0), extent=[(xmin-xshape)*pixelwidthPC,(xmax-xshape)*pixelwidthPC, \
		   (ymin-yshape)*pixelwidthPC,(ymax-yshape)*pixelwidthPC], origin='lower')
	ax1.set_xlabel('Distance from Centre in x-direction (pc)')
	ax1.set_ylabel('Distance from Centre in y-direction (pc)')
	ax1.set_title('T_max Map')
	### /Map

	#ax2 = fig.add_subplot(132)
	### Surface
	ax2.imshow(xi[0], interpolation = 'none', extent = [-dX*pixelwidthPC,dX*pixelwidthPC,-dY*pixelwidthPC,dY*pixelwidthPC], vmin=0, vmax=xi.max(), aspect='auto', origin='lower')
	levels = np.array([0.2,0.4,0.6,0.8])*xi[0].max()
	if xi[0].max() > 0:
		ax2.contour(xi[0], levels, extent=[-dX*pixelwidthPC,dX*pixelwidthPC,-dY*pixelwidthPC,dY*pixelwidthPC], vmin=0, vmax=xi.max(), colors='k')
	else:
		if xi[0].max() == 0:
			print xi[0].min()
	ax2.set_title('xi at 0 km/s')
	ax2.set_xlabel('Distance from Initial Location in x-direction (pc)')
	ax2.set_ylabel('Distance from Initial Location in y-direction (pc)')
	### /Surface

	#ax3 = fig.add_subplot(133)
	### Plot
	for i in range (0, np.abs(dV/ddV)+1):
	    if dV*velocityres > 0:
		ax3.plot(X, corr_funct[i],label='xi at +'+str('{0:.2f}'.format(i*ddV*velocityres*dV/np.abs(dV)))+' km/s')
	    elif dV*velocityres < 0:
		ax3.plot(X, corr_funct[i],label='xi at '+str('{0:.2f}'.format(i*ddV*velocityres*dV/np.abs(dV)))+' km/s')
	    else:
		ax3.plot(X, corr_funct[i],label='xi at '+str('{0:.2f}'.format(i*ddV*velocityres))+' km/s')
	ax3.set_title('Avg. Corr. Funct. vs. Radial "Distance" from Center of xi Plots')
	ax3.set_xlabel('Distance from Initial Location (pc)')
	ax3.set_ylabel('Average xi')
	ax3.legend(loc='upper right')
	### /Plot
	plt.tight_layout()
	plt.savefig('entire_xi_'+imagename+'.png')
	plt.clf()
