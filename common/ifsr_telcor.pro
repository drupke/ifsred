; docformat = 'rst'
;
;+
;
; Rebin data cube.
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Rebinned data cube.
;
; :Params:
;    infile: in, required, type=string
;      Path and filename of input file.
;    outfile: in, required, type=string
;      Path and filename of output file.
;    gcfile: in, required, type=string
;      Path and filename of telluric correction spectrum.
;
; :Keywords:
;    amcor: in, optional, type=double
;    wvcor: 
;    varnorm: in, optional, type=double
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
;      2010jun01, DSNR, created
;      2011may27, DSNR, fixed variance treatment; added AM correction
;      2014feb07, DSNR, complete rewrite for ease of use/customization;
;                       added detailed documentation
;
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
pro ifsr_telcor,infile,outfile,tcfile,amcor=amcor,varnorm=varnorm

  header=1
  cube = ifsf_readcube(infile,header=header,/quiet)  

  tcspec = readspec(tcfile)
  tcwave = tcspec[*,0]
  tcflux = tcspec[*,1]

  tcfluxi = interpol(tcflux,tcwave,cube.wave)
  
; To apply airmass correction, assume I = I_o exp(-amcor*tau).  Larger
; amcor means airmass in data is higher than in reference data, so
; correction must be larger.
  if keyword_set(amcor) then begin
     icor = where(tcfluxi lt 1)
     tcfluxi[icor] = alog(tcfluxi[icor])
     tcfluxi[icor] *= amcor
     tcfluxi[icor] = exp(tcfluxi[icor])
  endif

; Increase variance in telluric-corrected region to account for
; systematic uncertainty

  if ~ keyword_set(varnorm) then varnorm = 4d
; New = good
  tcfluxi_varnorm = tcfluxi
  incor = where(tcfluxi ge 1)
  icor = where(tcfluxi lt 1)
  tcfluxi_varnorm[incor] = 1d
  tcfluxi_varnorm[icor] = varnorm

  for i=0,cube.ncols-1 do begin
     for j=0,cube.nrows-1 do begin
        cube.dat[i,j,*] /= tcfluxi
        cube.var[i,j,*]  /= tcfluxi
        cube.var[i,j,*]  /= tcfluxi
        cube.var[i,j,*] *= tcfluxi_varnorm
     endfor
  endfor

; Write output

  writefits,outfile,cube.phu,header.phu
  writefits,outfile,cube.dat,header.dat,/append
  writefits,outfile,cube.var,header.var,/append
  writefits,outfile,cube.dq,header.dq,/append

end
