
# 5.13.16 - Calculates S_2 with noise correction using functions from "Cubes.py".

print('\nWelcome to Cubes_nonoise_multiple! \n \nAvailable functions: \n  image_make: Produces a 2D map and 1D plot of the noise-corrected S_2.  \n  multiple_images: Activates image_make for several preset subcubes all at once.\n \nThis program makes use of Cubes.py.\n \n')
import Cubes

def image_make(vmin, vmax, ymin, ymax, xmin, xmax, deltaX = 100, deltaV = 3, deltadeltaX = 10, filename="paws_norot", drawmap="False"):
	"""
	   Generates a noise-corrected 2D map and 1D plot of S_2, from subcube of the 
	   specified dimensions; using the .fits file in "Cubes.py".

	   Argument format: "(vmin,vmax, ymin,ymax, xmin,xmax, deltaX=100, deltaV=3,
	      deltadeltaX=10, filename="paws_norot", drawmap="False")."
	   ^ These are the parameters of the desired subcube, along with maximum dX/dY,
	     maximum dZ, "step size" for calculating S_2, the selected .fits file
	     name (minus the ".fits" extension), and an option to draw the subcube
	     as it appears on the galaxy T_max map (within the range specified).

	   Note that vmin and vmax will only affect the overall structure function (from
	     signal+noise), but not the noise-only structure function.


	   WARNING: Selecting too large of a subcube will hugely increase processing time.
	   If you use a large cube, be sure to set deltadeltaX to be larger in structgen.
	   For reliability, make sure that deltadeltaX is a factor of deltaX.

	   Be aware that processing time will increase with large deltaX and deltaV 
	   values, but can dramatically decrease with larger deltadeltaX at the cost of
	   plot resolution."""

	imagename = "M51_"+str(vmin)+"to"+str(vmax)+"_"+str(ymin)+"to"+str(ymax)+"_"+str(xmin)+"to"+str(xmax)

	subcube = Cubes.cubegen(vmin,vmax,ymin,ymax,xmin,xmax,filename,drawmap, imagename)	# Will draw a map of the subcube if drawmap="True".
	noisecube = Cubes.cubegen(100,120,ymin,ymax,xmin,xmax,filename,"False", imagename)	# We don't want a map of the noise, since it's just... noise. Nothing to see there.

	S2 = Cubes.structgen(subcube,deltaX,deltaV,deltadeltaX,False)
	S2n = Cubes.structgen(noisecube,deltaX,deltaV,deltadeltaX,False)
	
	S_2 = S2 - S2n		# This is the structure function from Signal ONLY.

	Cubes.mapgen(S_2, deltaX, deltaV, imagename)
	Cubes.plotgen(S_2, deltaX, deltaV, deltadeltaX, imagename)

def multiple_images(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=10):
	"""Activates image_make with each of the .py file's subcube selections,
	   all under spectral range (vmin,vmax) with maximum dX/dY, maximum dV,
	   and "step size". Be aware that this function will forcibly draw
	   all maps and plots involved.

	   Argument format: "(vmin=40,vmax=80, deltaX=40, deltaV=3, deltadeltaX=5).

	   WARNING: Selecting too large of a vmax-vmin will hugely increase
	   processing time."""

	filename = "paws_norot"
	drawmap = "True"

	image_make(vmin,vmax,350,550,500,700,deltaX,deltaV,deltadeltaX,filename,drawmap)		# Chunk of one of M51's arms, away from center.
#	image_make(vmin,vmax,200,400,425,625,deltaX,deltaV,deltadeltaX,filename,drawmap)		# Galaxy arm near center.
#	image_make(vmin,vmax,220,420,260,460,deltaX,deltaV,deltadeltaX,filename,drawmap)		# Galaxy arm near center (other side).
#	image_make(vmin,vmax,350,550,120,320,deltaX,deltaV,deltadeltaX,filename,drawmap)		# Chunk of galaxy arm, away from center.
#	image_make(vmin,vmax,350,550,250,450,deltaX,deltaV,deltadeltaX,filename,drawmap)		# Chunk of galaxy arm, not near center but not quite on outskirts either.
#	image_make(vmin,vmax,100,300,570,770,deltaX,deltaV,deltadeltaX,filename,drawmap)		# Big "empty" area.
#	image_make(vmin,vmax,200,400,360,560,deltaX,deltaV,deltadeltaX,filename,drawmap)		# Center of galaxy. (mostly centered on the nucleus)
