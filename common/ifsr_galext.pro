; docformat = 'rst'
;
;+
;
; Correct for Galactic extinction using the Fitzpatrick et al. 1999 curve.
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Corrected spectrum or spectra.
;
; :Params:
;    infile: in, required, type=string
;      Path and filename of input file.
;    outfile: in, required, type=string
;      Path and filename of output file.
;    ebv: in, required, type=double
;      
; :Keywords:
;    datext: in, optional, type=integer, default=1
;      Extension # of data plane. Set to a negative number if the correct
;      extension is 0, since an extension of 0 ignores the keyword.
;    varext: in, optional, type=integer, default=2
;      Extension # of variance plane.
;    dqext: in, optional, type=integer, default=3
;      Extension # of DQ plane. Set to a negative number if there is no DQ.
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
;      2022jan05, DSNR, created
;
; :Copyright:
;    Copyright (C) 2022 David S. N. Rupke
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
pro ifsr_galext,infile,outfile,ebv,datext=datext,varext=varext,dqext=dqext

   if ~keyword_set(datext) then datext=1
   if ~keyword_set(varext) then varext=2
   if ~keyword_set(dqext) then dqext=3
   
   header=1
   cube = ifsf_readcube(infile,header=header,datext=datext,varext=varext,$
      dqext=dqext)

   fm_unred,cube.wave,dblarr(cube.nz)+1d,ebv,fluxcor
   fluxcor = double(fluxcor)
   if cube.ndim eq 3 then $
      fluxcorarr = rebin(reform(fluxcor,1,1,cube.nz),cube.ncols,cube.nrows,cube.nz) $
   else if cube.ndim eq 2 then $
      fluxcorarr = rebin(reform(fluxcor,1,cube.nz),cube.ncols,cube.nz) $
   else $
      fluxcorarr = fluxcor

   newdat = cube.dat*fluxcorarr
   newvar = cube.var*fluxcorarr*fluxcorarr

   ;  Write output
   if datext gt 0 then begin
      writefits,outfile,cube.phu,header.phu
      writefits,outfile,newdat,header.dat,/append
   endif else begin
      writefits,outfile,newdat,header.dat
   endelse
   writefits,outfile,newvar,header.var,/append
   if dqext ne -1 then $
      writefits,outfile,cube.dq,header.dq,/append

end
