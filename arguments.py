import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as ss
import math

def summary(Array,ignore_nan=False):
	
	if ignore_nan==True:
		Mean = np.nanmean(Array)
		STD = np.nanstd(Array)
		Max = np.nanmax(Array)
	else:
		Mean = np.mean(Array)
		STD = np.std(Array)
		Max = np.max(Array)		

	return(Mean,STD,Max)
