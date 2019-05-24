; docformat = 'rst'
;
;+
;
; Read DR14 version of SDSS spectrum and write linearized spectrum into
; usual data, variance, DQ multi-extension FITS file. Extension 3 also includes
; spectral resolution in units of R = lam/dlam.
;
; :Categories:
;    IFSRED/SDSS
;
; :Returns:
;
; :Params:
;
; :Keywords:
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
;      2019mar27  DSNR  created
;
; :Copyright:
;    Copyright (C) 2019 David S. N. Rupke
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
;    along with this program.  If not, see http://www.gnu.org/licenses/.
;
;-
pro ifsr_ingestsdss,infile,outfile

   dummy = mrdfits(infile,0,header,/silent)
   instr = mrdfits(infile,1,/silent)
   
   dat = instr.flux
   var = 1d/instr.ivar
   wave = 10d^instr.loglam
   dq = instr.flux*0d
   nz = n_elements(dat)
   
;  According to
; https://data.sdss.org/datamodel/files/BOSS_SPECTRO_REDUX/RUN2D/spectra/PLATE4/spec.html
;  wdisp = Wavelength dispersion (sigma of fitted Gaussian) in units of number of pixel
   wdisp = instr.wdisp
;  convert to R = lam/dlam
   sres = wave[1:nz-1]/(wdisp*(wave[1:nz-1] - wave[0:nz-2])*2.35d)
   sres = [sres[0],sres]

;  Convert to air wavelengths
; https://www.sdss.org/dr14/spectro/spectro_basics/
   wave = wave / ( 1d +  5.792105d-2/(238.0185d - (1d4/wave)^2) + 1.67917d-3/( 57.362d - (1d4/wave)^2))
   
   waveold = wave
   datold = dat
   varold = var
   dqold = dq
   crpix = 1d
   cdelt = double(wave[nz-1]-wave[0])/double(nz-1)
   wave = double(wave[0]) + dindgen(nz)*cdelt
   crval = wave[0]
   dat = interpol(datold,waveold,wave,/spline)
   var = interpol(varold,waveold,wave,/spline)
   dq = interpol(dqold,waveold,wave)
;         print,'IFSF_READCUBE: Interpolating DQ; values > 0.01 set to 1.'
   ibd = where(dq gt 0.01,ctbd)
   if ctbd gt 0 then dq[ibd] = 1

    newheader = header
    sxaddpar,newheader,'CRVAL1',crval
    sxaddpar,newheader,'CRPIX1',crpix
    sxaddpar,newheader,'CDELT1',cdelt
    sxaddpar,newheader,'C1_1',cdelt

    mwrfits,dat,outfile,newheader,/create
    mwrfits,var,outfile,newheader,/silent
    mwrfits,dq,outfile,newheader,/silent
    mwrfits,sres,outfile,newheader,/silent

end
