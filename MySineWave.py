import matplotlib.pyplot as plt
import numpy as np
import sys

def MySineWave(Lam):

	L = Lam		# This is the wavelength.

	x = np.linspace(-2.0*np.pi, 2.0*np.pi, 200)
	y = np.sin(2.0*np.pi*x/L)

	plt.plot(x/np.pi,y,'k:')

	plt.title("Sine Wave")
	plt.xlabel("x (in units of pi)")
	plt.ylabel("y = sin(2*pi/L * x)")
	plt.axis([-2, 2, -1, 1])
	plt.show()

