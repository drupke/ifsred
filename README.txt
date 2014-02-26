-------------------------------------------------------------------------
OVERVIEW
-------------------------------------------------------------------------

IFSRED is a general-purpose library for reducing data from integral
field spectrographs (IFSs). For a general IFS data cube, it contains
routines to do the following:

- find and apply a zero-point shift in a wavelength solution on a
  spaxel-by-spaxel basis, using sky lines [in conjunction with IFSFIT
  -- see requirements below]
- find the spatial coordinates of a flux peak
- empirically correct for differential atmospheric refraction
- mosaic dithered exposures 
- (integer) rebin
- apply a telluric correction

IFSRED also contains a sky-subtraction routine designed specifically
to deal with GMOS data; however, this can be easily modified for any
instrument.

Finally, IFSRED contains some GMOS-specific software:

- Drop-in replacements for GFEXTRACT, GFTRANSFORM, GMOSAIC, and
  GSCALIBRATE in the v1.12 release of the Gemini IRAF package. These
  permit the propagation of a variance plane.
- A routine to create a map of the GMOS image plane prior to
  resampling from hexagonal to square pixels.

-------------------------------------------------------------------------
REQUIREMENTS
-------------------------------------------------------------------------

IDL (tested on v8.3, but probably backwards compatible with some older
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

-------------------------------------------------------------------------
USAGE
-------------------------------------------------------------------------

The IRAF routines in this package (in the subdirectory 'gmos/iraf/')
can be used by simply copying them to the proper directory. For a
Ureka installation, this directory is

/Your/path/Ureka/variants/common/iraf/gemini/gmos/

For help with the individual IDL routines, see the documentation of
each routine.

The EXAMPLE.txt file walks through the data reduction of a full data
cube, using Gemini IRAF tasks, the routines included in this package,
and the routines included in the IFSFIT package (see requirements).

-------------------------------------------------------------------------
QUESTIONS? BUGS? WANT TO MODIFY THE CODE?
-------------------------------------------------------------------------

Feel free to contact David Rupke at drupke@gmail.com with questions,
bug reports, etc.

Modifications are encouraged, but subject to the license.

-------------------------------------------------------------------------
LICENSE AND COPYRIGHT
-------------------------------------------------------------------------

Copyright (C) 2014 David S. N. Rupke

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
