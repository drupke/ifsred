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
;    tcfile: in, required, type=string
;      Path and filename of telluric correction spectrum.
;
; :Keywords:
;    amcor: in, optional, type=double, default=1
;      Multiplier to the absorption optical depth to either increase
;      or decrease correction. Higher amcor means airmass in data is
;      higher than in the reference data. This is equal to the
;      airmass ratio of the data and reference data.
;    varnorm: in, optional, type=double, default=1
;      Multiplier to the output variance to bias subsequent fits away
;      from telluric regions. The multiplier is applied to regions
;      where the input telluric correction spectrum is <1.
;    wvcor: in, optional, type=double, default=1
;      Multiplier to the water vapor absorption optical depth in regions of water
;      vapor absorption to account for changes in water vapor optical depth.
;      This is applied in addition to AMCOR if both are specified, and is applied
;      in the same way.
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
;      2014may01, DSNR, updated documentation
;      2015jan21, DSNR, added option to adjust water vapor corrections
;      2015jan23, DSNR, changed default varnorm to 1 (was 4)
;      2015may20, DSNR, re-wrote logic for dealing with cubes -- no longer
;                       assume 3-d data.
;
; :Copyright:
;    Copyright (C) 2014-2015 David S. N. Rupke
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
pro ifsr_telcor,infile,outfile,tcfile,amcor=amcor,varnorm=varnorm,wvcor=wvcor

  header=1
  cube = ifsf_readcube(infile,header=header,/quiet)  

  tcspec = ifsr_readspec(tcfile)
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

; Correction to decrease how corrections are applied in regions affected by 
; water vapor. Correction is applied in the same way as AMCOR.
  if keyword_set(wvcor) then begin
     icor = where(tcfluxi lt 1 AND $
                  ((tcwave gt 7100d AND tcwave lt 7500) OR $
                   (tcwave gt 7800d AND tcwave lt 8300)))
     tcfluxi[icor] = alog(tcfluxi[icor])
     tcfluxi[icor] *= wvcor
     tcfluxi[icor] = exp(tcfluxi[icor])
  endif

; Increase variance in telluric-corrected region to account for
; systematic uncertainty

  if ~ keyword_set(varnorm) then varnorm = 1d
; New = good
  tcfluxi_varnorm = tcfluxi
  incor = where(tcfluxi ge 1)
  icor = where(tcfluxi lt 1)
  tcfluxi_varnorm[incor] = 1d
  tcfluxi_varnorm[icor] = varnorm

  sizecube = size(cube.dat)

  if cube.ncols eq 1 then begin
     tcfluxi_use = tcfluxi
     tcfluxi_varnorm_use = tcfluxi_varnorm
  endif else if cube.ncols gt 1 then begin
     if cube.nrows gt 1 then begin
        tcfluxi_use = rebin(reform(tcfluxi,1,1,cube.nz),$
                            cube.ncols,cube.nrows,cube.nz)
        tcfluxi_varnorm_use = rebin(reform(tcfluxi_varnorm,1,1,cube.nz),$
                                    cube.ncols,cube.nrows,cube.nz)
     endif else begin
        if wavedim eq 1 then begin
           tcfluxi_use = rebin(reform(tcfluxi,cube.nz,1),$
                               cube.nz,cube.ncols)
           tcfluxi_varnorm_use = rebin(reform(tcfluxi_varnorm,cube.nz,1),$
                                       cube.nz,cube.ncols)
        endif else begin         
           tcfluxi_use = rebin(reform(tcfluxi,1,cube.nz),$
                               cube.ncols,cube.nz)
           tcfluxi_varnorm_use = rebin(reform(tcfluxi_varnorm,1,cube.nz),$
                                       cube.ncols,cube.nz)
        endelse
     endelse
  endif
  
  cube.dat /= tcfluxi_use
  cube.var /= tcfluxi_use
  cube.var /= tcfluxi_use
  cube.var *= tcfluxi_varnorm_use

; Write output

  writefits,outfile,cube.phu,header.phu
  writefits,outfile,cube.dat,header.dat,/append
  writefits,outfile,cube.var,header.var,/append
  writefits,outfile,cube.dq,header.dq,/append

end
