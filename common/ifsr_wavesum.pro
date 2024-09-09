; docformat = 'rst'
;
;+
;
; Sum data cube over a wavelength range.
;
; :Categories:
;    IFSRED
;
; :Returns:
;    X x Y x 2 array (with X and Y spatial dimensions) containing data
;    and variance.
;
; :Params:
;    cube: in, required, type=structure
;      Structure with input data cube; output of IFSF_READCUBE.
;    wavesum: in, required, type=dblarr(2)
;      Lower and upper limits of wavelength region over which to sum.
;
; :Keywords:
;    applydq: in, optional, type=bool
;      Apply DQ plane.
;    average: in, optional, type=bool
;      Calculate average flux/A instead of total flux. If there are
;      missing pixels, the total flux calculation adds back in the
;      missing pixels as having the average flux.
;    npix: out, optional, type=int
;      Maximum number of pixels summed over wavelength.
;
; :Author:
;    David S. N. Rupke::
;      Rhodes College
;      Department of Physics
;      2000 N. Parkway
;      Memphis, TN 38104
;      drupke@gmail.com
;
; :History:
;    ChangeHistory::
;      2018nov04, DSNR, created
;
; :Copyright:
;    Copyright (C) 2018 David S. N. Rupke
;
;    This program is free software: you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation, either version 3 of
;    the License or any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program.  If not, see
;    http://www.gnu.org/licenses/.
;
;-
function ifsr_wavesum,cube,wavesum,average=average,npix=npix,applydq=applydq

   bad=1d99

;  keep state of input cube
   incube = cube

;  Locations of wavelengths to sum between
   iwavesum = value_locate(cube.wave,wavesum)
;  Total # of pixels summed in spectral dimension
   maxwavepix = iwavesum[1]-iwavesum[0]+1
   totnwavepix = intarr(cube.ncols,cube.nrows)+maxwavepix
;  Take care of pixels with DQ > 0. Assume DQ is either 0 (good) or 1 (bad).
   if keyword_set(applydq) then begin
;     Total value of DQ over spectral dim.
      totdq = total(cube.dq[*,*,iwavesum[0]:iwavesum[1]],3)
;     spaxels where total DQ > 0
      ibaddq = where(totdq gt 0,ctbaddq)
      if ctbaddq gt 0 then begin
;        Any pixel (spatial or spectral) where DQ > 0
         ibaddqcube = where(cube.dq gt 0)
;        Set data = 0 where DQ > 0
         cube.dat[ibaddqcube] = 0d
;        Subtract total DQ = 1 pixels in spectral dim. from the total # of
;        spectral pixels that contribute to the sum in each spaxel
         totnwavepix[ibaddq] -= totdq[ibaddq]
      endif
   endif
;  Sum data in spectral dim. Note that this does not interpolate over DQ > 0 pixels,
;  just ignores them.
   totflux = total(cube.dat[*,*,iwavesum[0]:iwavesum[1]],3,/double)
   totvar = total(cube.var[*,*,iwavesum[0]:iwavesum[1]],3,/double)
;  Set points with no data to 0 flux and bad variance
   outflux = dblarr(cube.ncols,cube.nrows)
   outvar = dblarr(cube.ncols,cube.nrows)+bad
   igddat = where(totnwavepix gt 0,ctgddat)
;  Calculate average flux density per spectral pixel
   if keyword_set(average) then begin
      outflux[igddat] = totflux[igddat] / double(totnwavepix[igddat])
      outvar[igddat] = totvar[igddat] / double(totnwavepix[igddat])
;  ... or make a correction for missing pixels. Note that this is *not* an 
;  interpolation.
   endif else begin
;     Multiply by dispersion to convert from flux density to flux
      outflux[igddat] *= cube.cdelt
      outvar[igddat] *= cube.cdelt^2d
      outflux[igddat] *= double(maxwavepix)/double(totnwavepix[igddat])
      outvar[igddat] *= double(maxwavepix)/double(totnwavepix[igddat])
   endelse

;  Restore state of input cube
   cube = incube
   
;  Return max. number of pixels in average
   npix = maxwavepix

   return,[[[outflux]],[[outvar]]]

end
