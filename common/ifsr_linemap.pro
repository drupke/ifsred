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
;      2018nov04, DSNR, created
;      2019jan25, DSNR, separated IFSR_MAKELINEMAP into subroutine
;
; :Copyright:
;    Copyright (C) 2018-19 David S. N. Rupke
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
function ifsr_makelinemap,cube,wavesum,wavesub=wavesub

   ;  Get linemap
   linesum = ifsr_wavesum(cube,wavesum,/average,npix=npix,/applydq)
   linesumdat = linesum[*,*,0]
   linesumvar = linesum[*,*,1]
   igdline = where(linesumdat gt 0,ctgdline)
   ;  Subtract continuum if requested
   if keyword_set(wavesub) then begin
      subsumdat = dblarr(cube.ncols,cube.nrows)
      subsumvar = dblarr(cube.ncols,cube.nrows)
      subsumdq = dblarr(cube.ncols,cube.nrows)
      sizews = size(wavesub)
      if sizews[0] eq 1 then nws=1
      if sizews[0] eq 2 then nws=sizews[2]
      for i=0,nws-1 do begin
         subsum_tmp = ifsr_wavesum(cube,wavesub[*,i],/average,/applydq)
         subsumdat_tmp = subsum_tmp[*,*,0]
         subsumvar_tmp = subsum_tmp[*,*,1]
         igddat = where(subsumdat_tmp gt 0,ctgddat)
         subsumdat[igddat] += subsumdat_tmp[igddat]
         subsumvar[igddat] += subsumvar_tmp[igddat]
         subsumdq[igddat] += 1d
      endfor
      ;     average over number of regions
      igddat = where(subsumdq gt 0d,ctgddat)
      subsumdat[igddat] /= subsumdq[igddat]
      subsumvar[igddat] /= subsumdq[igddat]
      ;     subtract
      linesumdat[igdline] -= subsumdat[igdline]
      linesumvar[igdline] += subsumvar[igdline]
   endif
   ;  ... and return to total flux rather than average per A. Also
   ;  multiply by dispersion to convert from flux/A to flux
   linesumdat[igdline] *= double(npix)*cube.cdelt
   linesumvar[igdline] *= double(npix)*cube.cdelt^2d

   linesum = {dat: linesumdat, var: linesumvar}
   
   return,linesum

end

pro ifsr_linemap,infile,outfile,wavesum,wavesub=wavesub,$
                 datext=datext,varext=varext,dqext=dqext

   bad=1d99

   if ~ keyword_set(datext) then datext=1
   if ~ keyword_set(varext) then varext=2
   if ~ keyword_set(dqext) then dqext=3

   header=1
   cube = ifsf_readcube(infile,header=header,datext=datext,varext=varext,$
                        dqext=dqext,/quiet)

   linesum = ifsr_makelinemap(cube,wavesum,wavesub=wavesub)

   appenddat=0b
   if datext eq 1 then begin
      writefits,outfile,[],header.phu
      appenddat=1b
   endif
   writefits,outfile,linesum.dat,header.dat,append=appenddat
   writefits,outfile,linesum.var,header.var,/append

end
