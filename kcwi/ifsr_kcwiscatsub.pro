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
;      2018aug11  DSNR  created
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
;    along with this program.  If not, see http://www.gnu.org/licenses/.
;
;-
pro ifsr_kcwiscatsub,instr,gapfile,scatfun,parinfo,$
                     dy=dy,psfcent=psfcent,$
                     plots=plots,argscatfun=argscatfun,$
                     sigclip=sigclip

   if ~ keyword_set(dy) then dy=100
   if ~ keyword_set(plots) then plots=0 else plots=1

;  Gap file
   readcol,gapfile,x0,x1,format='(I0,I0)',/silent
   igdx = []
   for j=0,n_elements(x0)-1 do $
      igdx = [igdx,indgen(x1[j]-x0[j]+1)+x0[j]-1]

;  In
   int = mrdfits(instr+'_int.fits',0,ihead,/silent)
   var = mrdfits(instr+'_var.fits',0,vhead,/silent)

;  Number of bins in y direction
   isize = size(int)
   nbin = floor(double(isize[2])/double(dy))+1

;  Bad x-values
   iallx = indgen(isize[1])
   ibdx = cgsetdifference(iallx,igdx)

;  x array for fitting, in single-offset values
   xbin = dindgen(isize[1])+1

;  y arrays for fitting
   ybin = dblarr(isize[1])
   ybinvar = dblarr(isize[1])
   ybinerr = dblarr(isize[1])
;  y fits to all bins
   ysurf = dblarr(isize[1],nbin)
;  central coordinate of bin in single-offset coordinates
   ybinloc = dblarr(nbin)

   if keyword_set(plots) then set_plot,'x'
;  Loop through y bins
   for i=0,nbin-1 do begin
;     Sum data
      iystart=i*dy
      if i ne nbin-1 then iyend = (i+1)*dy-1 else iyend = isize[2]-1
;      ybin=total(int[*,iystart:iyend],2,/double)/(iyend-iystart+1)
      ybin=median(int[*,iystart:iyend],dim=2,/double)
      ybinloc[i]=double(iyend-iystart)/2d + double(iystart)
      ybinvar=total(var[*,iystart:iyend],2,/double)/(iyend-iystart+1)
      ybinerr=sqrt(ybinvar)/sqrt(iyend-iystart+1) ; std err
;     Set values not to fit to NaN
      ybin[ibdx] = !VALUES.D_NAN
      ybinerr[ibdx] = !VALUES.D_NAN
      if keyword_set(argscatfun) then $
         coeffs = mpfitfun(scatfun,xbin,ybin,ybinerr,/nan,yfit=yfit,$
                           parinfo=parinfo,functargs=argscatfun,/quiet) $
      else $
         coeffs = mpfitfun(scatfun,xbin,ybin,ybinerr,/nan,yfit=yfit,$
                           parinfo=parinfo,/quiet)
      if plots then begin
         cgplot,xbin,ybin,psym=16,err_ylo=ybinerr,err_yhi=ybinerr,/err_clip,$
                err_col='Red'
         cgoplot,xbin,yfit
      endif
;      3-sigma clip
      ctclip=0
      if keyword_set(sigclip) then begin
         ydiff = abs(yfit-ybin)/ybinerr
         iclip = where(ydiff ge sigclip,ctclip)
         if ctclip gt 0 then begin
            xclip = xbin[iclip]
            yclip = ybin[iclip]
            ybin[iclip] = !VALUES.D_NAN
            ybinerr[iclip] = !VALUES.D_NAN
            if keyword_set(argscatfun) then $
               coeffs = mpfitfun(scatfun,xbin,ybin,ybinerr,/nan,yfit=yfit,$
                                 parinfo=parinfo,functargs=argscatfun,/quiet) $
            else $
               coeffs = mpfitfun(scatfun,xbin,ybin,ybinerr,/nan,yfit=yfit,$
                                 parinfo=parinfo,/quiet)
         endif
      endif
      ysurf[*,i] = yfit
      if plots then begin
         if ctclip gt 0 then cgoplot,xclip,yclip,psym=7,symsize=2,color='Blue'
         q = ''
         read,'Next? (Q-quit prompting for plots, <cr>-next): ',q
         if strupcase(strmid(q,0,1)) eq 'Q' then plots = 0
       endif
   endfor

;  Now interpolate results in between
   ysurfit = dblarr(isize[1],isize[2])
   for i=0,isize[1]-1 do $
      ysurfit[i,*] = interpol(ysurf[i,*],ybinloc,dindgen(isize[2])+1,/spline)

;  Subtract scattered light
   int -= ysurfit

;  Out
   writefits,instr+'_intd.fits',int,ihead

end
