; docformat = 'rst'
;
;+
;
; Read single spectrum in FITS format.
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Wavelength and flux arrays of spectrum.
;
; :Params:
;    infile: in, required, type=string
;      Path and filename of input  FITS file.
;
; :Keywords:
;    extension: in, optional, type=integer
;      Which extension to read.
;    header: out, optional, type=strarr
;      Optionally, output header.
;    waveext: in, optional, type=integer
;      The extention number of a wavelength array.
;    
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
;      2009aprXY, DSNR, created (because there is nothing like it in Astrolib)
;      2009jul01, DSNR, modified to include log-linear coordinates
;      2011jan21, DSNR, added option to output header
;      2013mar05, DSNR, modified treatment of error messages; not sure how
;                     this will work ...
;      2016aug24, DSNR, ported to IFSRED
;      2018feb08, DSNR, added WAVEEXT keyword
;
;
; :Copyright:
;    Copyright (C) 2016--2018 David S. N. Rupke
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
function ifsr_readspec,infile,extension=extension,header=header,waveext=waveext

  errmsg=''
  
  if ~ keyword_set(extension) then extension=0

  flux = mrdfits(infile,extension,header,/silent)

  if ~ keyword_set(waveext) then begin
     s = size(flux)
     if s[0] eq 1 then naps=1 else naps=s[2]
     npix = n_elements(flux)/naps
     wave = dindgen(npix)
     ct = intarr(6)
     ctype = sxpar(header,'CTYPE1',/silent,count=count)
     ct[0] = count
     crval = double(sxpar(header,'CRVAL1',/silent,count=count))
     ct[1] = count
     crpix = double(sxpar(header,'CRPIX1',/silent,count=count))
     ct[2] = count
     cdelt = double(sxpar(header,'CDELT1',/silent,count=count))
     ct[3] = count
     cd11  = double(sxpar(header,'CD1_1',/silent,count=count))
     ct[4] = count
     dcflag= sxpar(header,'DC-FLAG',/silent,count=count)
     ct[5] = count

     disp = cdelt
     nopar = where(ct eq 0,ctpar)
     if ctpar then begin
        if ctpar gt 1 then begin
           errmsg = 'Header keywords missing.'
           goto,error
        endif else begin
           if nopar[0] eq 3 then disp = cd11
        endelse
     endif

     if (ctype eq 'LINEAR  ' OR ctype eq 'lambda  ') then begin
        wave = crval + disp*(wave-crpix+1)
        if dcflag then wave = 10d^wave
     endif else begin
        errmsg = 'Dispersion solution not linear.'
        goto,error
     endelse
  endif else begin
     wave = mrdfits(infile,waveext,header_wave,/silent)
  endelse

error:
  if errmsg ne '' then print,'IFSR_READSPEC: ERROR: ',errmsg,$
                             format='(T5,A,A)'

  return,[[wave],[flux]]

end
