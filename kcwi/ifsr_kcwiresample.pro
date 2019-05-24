; docformat = 'rst'
;
;+
;
; :Categories:
;    IFSRED/KCWI
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
;      2018aug10  DSNR  created
;      2019mar28  DSNR  fixed resampling scale to keep constant total flux
;
; :Copyright:
;    Copyright (C) 2018--2019 David S. N. Rupke
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
pro ifsr_kcwiresample,instr,outfile,goodreg

;  Raw data
   dat = mrdfits(instr+'_icubes.fits',0,ihead,/silent)
   var = mrdfits(instr+'_vcubes.fits',0,vhead,/silent)
   msk = double(mrdfits(instr+'_mcubes.fits',0,mhead,/silent))
;  Data sections
   pxscl = sxpar(ihead,'PXSCL',/silent)*3600d
   slscl = sxpar(ihead,'SLSCL',/silent)*3600d

   datsize = size(dat)
;  Sizes of x and y spaxel regions specified on commnand line
   dx = goodreg[2]-goodreg[0]+1
   dy = goodreg[3]-goodreg[1]+1
   dz = datsize[3]
   
;  Number of output x spaxels is current number times ratio of pixel sizes. 
;  We round to nearest pixel
   newdx = floor(slscl/pxscl*double(dx))+1
   newdat = dblarr(newdx,dy,dz)
   newvar = dblarr(newdx,dy,dz)
   newmsk = dblarr(newdx,dy,dz)
;  positions in old y array.
   yarr = dindgen(dx)+0.5
;  positions in new yarr
   newyarr = (dindgen(newdx)+0.5)*pxscl/slscl
   
   for i=0,dz-1 do begin
      for j=0,dy-1 do begin
         newdat[*,j,i] = $
            interpol(dat[goodreg[0]-1:goodreg[2]-1,goodreg[1]-1+j,i],$
                     yarr,newyarr)*pxscl/slscl
         newvar[*,j,i] = $
            interpol(var[goodreg[0]-1:goodreg[2]-1,goodreg[1]-1+j,i],$
                     yarr,newyarr)*pxscl/slscl
         newmsk[*,j,i] = $
            interpol(msk[goodreg[0]-1:goodreg[2]-1,goodreg[1]-1+j,i],$
                     yarr,newyarr)
      endfor
   endfor

   ibad = where(newmsk gt 0)
   newmsk[ibad] = 1
   newmsk = fix(newmsk)

;  Adjust CD matrix for resampling.
;  This only works for PA=0 right now!
;  To get it working for other PAs, have to interpret CD matrix properly.
;   ra = sxpar(ihead,'CRVAL1',/silent)
   pa = sxpar(ihead,'ROTPOSN',/silent)
   if fix(pa) eq 0 then begin
      print,'IFSR_KCWIRESAMPLE: PA = 0, so adjusting CD1_1 for resampling.'
      oldcdx = sxpar(ihead,'CD1_1',/silent)
      oldcdy = sxpar(ihead,'CD2_2',/silent)
      newcdy = oldcdy
      newcdx = oldcdx*double(dx)/double(newdx)
      sxaddpar,ihead,'CD1_1',newcdx
      sxaddpar,ihead,'CD2_2',newcdy
      sxaddpar,vhead,'CD1_1',newcdx
      sxaddpar,vhead,'CD2_2',newcdy
      sxaddpar,mhead,'CD1_1',newcdx
      sxaddpar,mhead,'CD2_2',newcdy
   endif

   mwrfits,newdat,outfile,ihead,/create
   mwrfits,newvar,outfile,vhead
   mwrfits,newmsk,outfile,mhead

end
