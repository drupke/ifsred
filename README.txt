-------------------------------------------------------------------------
OVERVIEW
-------------------------------------------------------------------------

IFSRED is a general-purpose library for reducing data from integral
field spectrographs (IFSs). For a general IFS data cube, it contains
routines to do the following:

- apply a telluric correction [IFSR_TELCOR]
- empirically correct for differential atmospheric refraction
- find and apply a zero-point shift in a wavelength solution on a
  spaxel-by-spaxel basis, using sky lines [IFSF_SKY; in conjunction
  with IFSFIT -- see requirements below]
- find the spatial coordinates of a flux peak [IFSR_PEAK]
- flip/rotate a data cube [IFSR_FLIPROT]
- Gaussian smooth to seeing-match data [IFSR_BLUR]
- mosaic dithered exposures [IFSR_MOSAIC]
- rebin [IFSR_REBIN]
- sum over wavelength to make a linemap [IFSR_LINEMAP,
  IFSR_MAKELINEMAP, IFSR_WAVESUM]
- sum over spatial apertures into a single spectrum [IFSR_SPAXSUM]
- sky subtract [IFSR_SKYSUB]
- Voronoi bin [IFSR_VORBIN, IFSR_VORONOI]

IFSRED also contains various routines for instrument-specific
reduction or analysis:
- GMOS
- KCWI
- MaNGA
- NIFS
- OSIRIS
- SDSS
- WiFeS [S7]

-------------------------------------------------------------------------
REQUIREMENTS
-------------------------------------------------------------------------

IDL (tested on v8.3-v8.7, but probably backwards compatible with some older
versions)

IDL libraries:
- IFSFIT
  https://github.com/drupke/ifsfit
- IDL Astronomy User's Library
  http://idlastro.gsfc.nasa.gov
- MPFIT
  http://www.physics.wisc.edu/~craigm/idl/idl.html
- Coyote
  http://www.idlcoyote.com/documents/programs.php#COYOTE_LIBRARY_DOWNLOAD
  [or from the subversion repository: https://code.google.com/p/idl-coyote/]
- VORONOI_2D_BINNING
  https://www-astro.physics.ox.ac.uk/~mxc/software/#binning
  [for Voronoi binning]

-------------------------------------------------------------------------
USAGE
-------------------------------------------------------------------------

For help with individual IDL routines, see the documentation of each
routine.

-------------------------------------------------------------------------
QUESTIONS? BUGS? WANT TO MODIFY THE CODE?
-------------------------------------------------------------------------

Feel free to contact David Rupke at drupke@gmail.com with questions,
bug reports, etc.

Modifications are encouraged, but subject to the license.

-------------------------------------------------------------------------
LICENSE AND COPYRIGHT
-------------------------------------------------------------------------

Copyright (C) 2014--2019 David S. N. Rupke

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License or any
later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses/.
