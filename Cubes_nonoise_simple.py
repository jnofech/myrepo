
# 5.12.16 - Calculates S_2 with noise correction using functions from "Cubes.py".

print('\nWelcome to Cubes_nonoise_simple! \n \nAvailable functions: \n  image_make: Produces a 2D map and 1D plot of the noise-corrected S_2. \n \nThis program makes use of Cubes.py. \n \n')
import Cubes

def image_make(vmin, vmax, ymin, ymax, xmin, xmax, deltaX = 100, deltaV = 3, deltadeltaX = 10, deltadeltaV = 1, filename="paws_norot"):
	"""Generates a noise-corrected 2D map and 1D plot of S_2, from subcube of the 
	   specified dimensions; using the .fits file in "Cubes.py".

	   Argument format: "(vmin,vmax, ymin,ymax, xmin,xmax, deltaX=100, deltaV=3,
	      deltadeltaX=10, deltadeltaV=1, filename="paws_norot)."
	   ^ These are the parameters of the desired subcube, along with maximum dX/dY,
	     maximum dV, "step sizes" for calculating S_2, and the selected .fits file
	     name (minus the ".fits" extension).
	   Note that vmin and vmax will only affect the overall structure function (from
	     signal+noise), but not the noise-only structure function.


	   WARNING: Selecting too large of a subcube will hugely increase processing time.
	   If you use a large cube, be sure to set deltadeltaX to be larger in structgen.

	   Be aware that processing time will increase with large deltaX and deltaV 
	   values, but can decrease with larger deltadeltaX at the cost of plot
	   resolution (or with larger deltadeltaZ at the cost of the number of 1D
	   plots)."""

	
	subcube = Cubes.cubegen(vmin,vmax,ymin,ymax,xmin,xmax)
	noisecube = Cubes.cubegen(0,20,ymin,ymax,xmin,xmax)

	S2 = Cubes.structgen(subcube,deltaX,deltaV,deltadeltaX,deltadeltaV,False)
	S2n = Cubes.structgen(noisecube,deltaX,deltaV,deltadeltaX,deltadeltaV,False)
	
	S_2 = S2 - S2n		# This is the structure function from Signal ONLY.

	Cubes.mapgen(S_2, deltaX, deltaV, deltadeltaV)
	Cubes.plotgen(S_2, deltaX, deltaV, deltadeltaX, deltadeltaV)
