#! /usr/bin/env python

"""Stack 2D FITS image along columns and plot result.
   Assumes data is in first extension."""
   

import sys
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
 
if len(sys.argv) != 2:
   raise ValueError('Command-line argument is KCWI .fits file')

#print(f'Opening {sys.argv[1]}')
hdul = fits.open(sys.argv[1])
hdul.info()
spec = hdul[0].data
print(f'Stacking along columns')
stack = np.sum(spec,axis=0)
hdul.close()

plt.plot(stack)
plt.show()
