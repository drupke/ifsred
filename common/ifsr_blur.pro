; docformat = 'rst'
;
;+
;
; :Categories:
;    IFSRED
;
; :Returns:
;    
; :Params:
;    infile: in, required, type=strarr
;      Filenames for data cubes of individual exposures.
;    outfile: in, required, type=string
;      Path and filename of output file.
;    fwhm_gaussian: in, required, type=double
;      FWHM of Gaussian with which to convolve, in pixels.
;
; :Keywords:
;    nophu: in, optional, type=byte
;      Added flag to indicate that the 0th extension contains the data, not
;      the PHU.
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
;       2016dec06, DSNR, created
;       
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
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
pro ifsr_blur,infile,outfile,fwhm_gaussian,nophu=nophu

;  Loading data

   header=1
   if keyword_set(nophu) then begin
      datext=-1
      varext=1
      dqext=2
   endif else begin
      datext=1
      varext=2
      dqext=3
   endelse
   cube = ifsf_readcube(infile,header=header,/quiet,$
                        datext=datext,varext=varext,dqext=dqext)

;  Initialize arrays for new cube.
   datnew = dblarr(cube.ncols,cube.nrows,cube.nz)
   varnew = dblarr(cube.ncols,cube.nrows,cube.nz)

;  convert FWHM to box_width
   box_width = round(fwhm_gaussian^2d / 2d + 1d)

   ibad = where(cube.dq eq 1b,ctbad)
   if ctbad gt 0 then cube.dat[ibad] = !VALUES.D_NAN
   for i=0,cube.nz-1 do begin
      datnew[*,*,i] = filter_image(cube.dat[*,*,i],/all,/iter,smooth=box_width)
      varnew[*,*,i] = filter_image(cube.var[*,*,i],/all,/iter,smooth=box_width)
   endfor

   appenddat=1
   if ~ keyword_set(nophu) then writefits,outfile,cube.phu,header.phu $
   else appenddat=0
   writefits,outfile,datnew,header.dat,append=appenddat
   writefits,outfile,varnew,header.var,/append
   writefits,outfile,cube.dq,header.dq,/append

   datnew=0
   varnew=0

end
