; docformat = 'rst'
;
;+
;
; Sum emission line fluxes over wavelength to make a spatial map.
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Structure containing linemap and variance.
;
; :Params:
;    cube: in, required, type=structure
;      Structure with input data cube; output of IFSF_READCUBE.
;    wavesum: in, required, type=dblarr(2)
;      Lower and upper limits of wavelength region over which to sum.
;
; :Keywords:
;    allowneg: in, optional, type=boolean
;      Default is to calculate continuum fluxes using positive values
;      only; this allows negative fluxes in the computation as well.
;    wavesub: in, optional, type=dblarr(2,N)
;      Lower and upper wavelength regions for subtraction continuum.
;      E.g., [5000,5100] or [[4800,4900],[5000,5100]].
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
;      2019mar04, DSNR, moved IFSR_MAKELINEMAP to separate file from IFSR_LINEMAP
;      2019may22, DSNR, keyword ALLOWNEG allows negative fluxes in computing 
;                       continuum region averages to subtract
;      2024jan08, DSNR, ALLOWNEG now also allows good line regions to go neg
;
; :Copyright:
;    Copyright (C) 2019--2024 David S. N. Rupke
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
function ifsr_makelinemap,cube,wavesum,wavesub=wavesub,allowneg=allowneg

   bad = 1d99
   ;  Get linemap
   linesum = ifsr_wavesum(cube,wavesum,/average,npix=npix,/applydq)
   linesumdat = linesum[*,*,0]
   linesumvar = linesum[*,*,1]
   if keyword_set(allowneg) then begin
      igdline = where(linesumdat ne bad,ctgdline)
   endif else begin
      igdline = where(linesumdat gt 0,ctgdline)
      print,'IFSR_MAKELINEMAP: Applying multipliers/subtractions to positive fluxes only.'
   endelse
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
         if keyword_set(allowneg) then begin
            igddat = where(subsumdat_tmp ne 0,ctgddat)
         endif else begin
            igddat = where(subsumdat_tmp gt 0,ctgddat)
            print,'IFSR_MAKELINEMAP: Computing continuum region averages from positive fluxes only.'
         endelse
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
   ;  ... and return to flux density rather than average flux density per pixel. Also
   ;  multiply by dispersion to convert from flux/A to flux
   linesumdat[igdline] *= double(npix)*cube.cdelt
   linesumvar[igdline] *= double(npix)*cube.cdelt^2d

   linesum = {dat: linesumdat, var: linesumvar}

   return,linesum

end
