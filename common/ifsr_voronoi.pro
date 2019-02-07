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
pro ifsr_voronoi,inmap,center,targetSN,threshSN,outxdr,$
                 datext=datext,varext=varext ;,_extra=_extra

   bad=1d99
   if ~ keyword_set(datext) then datext=1
   if datext eq -1 then datext=0
   if ~ keyword_set(varext) then varext=2

   header=1
   linedat = readfits(inmap,ext=datext)
   linevar = readfits(inmap,ext=varext)
   sizearr = size(linedat)
   dx = sizearr[1]
   dy = sizearr[2]

   linesn = linedat/sqrt(linevar)
   igddat = where(linedat gt 0 AND linesn gt threshSN,ctgddat)
   lineerr = dblarr(dx,dy)+bad
   lineerr[igddat] = sqrt(linevar[igddat])

   print,'IFSR_VORONOI: There are '+string(ctgddat,format='(I0)')+$
         ' pixels with S/N > ',string(threshSN,format='(D0.2)')

;  Output locations in single-offset pixel coordinates
   xarr = indgen(dx)+1 - center[0]
   xvals = rebin(xarr,dx,dy)
   yarr = indgen(dy)+1 - center[1]
   yvals = rebin(reform(yarr,1,dy),dx,dy)

;  Load a colortable and open a graphic window
   set_plot,'x'
   loadct, 5
   r = GET_SCREEN_SIZE()
   window, xsize=r[0]*0.6, ysize=r[1]*0.8

   voronoi_2d_binning, xvals[igddat], yvals[igddat], linedat[igddat], $
                       lineerr[igddat], targetSN, $
                       binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale,$
                       pixsize=1,/plot

   vordat = {center: center,$
             binNum: binNum,$
             xin: xvals[igddat],$
             yin: yvals[igddat],$
             xNode: xNode,$
             yNode: yNode,$
             xBar: xBar,$
             yBar: yBar,$
             sn: sn,$
             nPixels: nPixels,$
             scale: scale}
   
   save,vordat,file=outxdr

end
