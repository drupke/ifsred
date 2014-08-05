; docformat = 'rst'
;
;+ 
; 
; Remove scattered light from NIFS data. Looks for minimum flux in gaps between
; spectrum blocks and fits a Moffat profile to row numbers of flux minima. Repeats
; in column blocks of width dx.
;
; :Categories:
;    IFSRED/NIFS
;
; :Returns:
;    Fits file with scattered light removed.
;
; :Params:
;    infile: in, required, type=string
;      Input MEF.
;    outfile: in, required, type=string
;      Output MEF.
;    ygap1: in, required, type=string
;      Row number of first spectrum.
;
; :Keywords:
;    plot: in, optional, type=boolean
;      Show plots of Moffat fits to flux minima in gaps between spectrum blocks.
;    lowlight: in, optional, type=boolean
;      Add extra constraints when there is little scattered light.
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
;      2013apr05, DSNR, created
;      2013apr23, DSNR, edited so in case of fitting scattered light below
;                       flux = 0, would not add back in flux; also added
;                       ygap1 as a parameter b/c of shifts with time/
;                       wavelength; and added "lowlight" case to better fit 
;                       when not much scattered light
;      2014aug05, DSNR, added documentation, license, and copyright
;
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
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
pro nifs_scatter,infile,outfile,ygap1,plot=plot,lowlight=lowlight
  
  if ~ keyword_set(plot) then doplot = 0 else doplot = 1

  sx = 2040 ; columns in input file
  sy = 2080 ; rows in input file
; column bin size for computing / removing scattered light
  if lowlight then dx = 100 else dx = 50
  nx = sx/dx ; number of bins

  dgap = 71d ; distance between gaps in spectral blocks
  ngap = 29 ; number of gaps between spectral blocks

; read data
  phu = mrdfits(infile,0,hdr_phu,ext=0,/silent)
  sci = mrdfits(infile,1,hdr_sci,ext=1,/silent)
  var = mrdfits(infile,2,hdr_var,ext=2,/silent)
  dq = mrdfits(infile,3,hdr_dq,ext=3,/silent)
  mdf = mrdfits(infile,4,hdr_mdf,ext=4,/silent)

  if doplot then set_plot,'x'

  ygap = dblarr(ngap) ; row number of minimum flux in gap
  fgap = dblarr(ngap) ; flux of minimum flux in gap
  ysec = dindgen(sy) ; row number

; In case of less scattered light, constrain profile first using
; entire spatial profile
  if lowlight then begin
     i1 = 0
     i2 = sx-1
     xsec = median(sci[i1:i2,*],dim=1)
     for j=0,ngap-1 do begin
        ygapcent = ygap1 + j*dgap
        iy1 = ygapcent-5
        iy2 = ygapcent+5
        xsecgap = xsec[iy1:iy2]
        fgap[j] = min(xsecgap,imingap)
        ygap[j] = iy1+imingap
     endfor
     fit = mpfitpeak(ygap,fgap,refparam,/moffat,nterms=6)
     ineg = where(fit lt 0,ctneg)
     if ctneg gt 0 then fit[ineg] = 0
     if doplot then begin
        cgplot,ygap,fgap,xran=[0,sx],yran=[min(fgap),max(fgap)],/xsty,/ysty
        cgplot,ygap,fit,/overplot,color='Blue'
        cgplot,ysec,xsec,/overplot,color='Red'
        wait,2
     endif
  endif

; Now cycle through chunks of image of width dx
  for i=0,nx do begin
     i1 = i*dx
     i2 = (i+1)*dx-1
     if i2 gt sx-1 then i2=sx-1
     xsec = median(sci[i1:i2,*],dim=1)

     for j=0,ngap-1 do begin
        ygapcent = ygap1 + j*dgap
        iy1 = ygapcent-5
        iy2 = ygapcent+5
        xsecgap = xsec[iy1:iy2]
        fgap[j] = min(xsecgap,imingap)
        ygap[j] = iy1+imingap
     endfor

;    Fix only the HWHM of the moffat fit.
     if keyword_set(lowlight) then begin
        parinfo = REPLICATE({$
                  limited:[0b,0b],$
                  limits:[0d,0d],$
                  value:0d,$
                  fixed:0b},6)
        parinfo.value = refparam
        parinfo[2].fixed=1b
        fit = mpfitpeak(ygap,fgap,param,/moffat,nterms=6,parinfo=parinfo)
     endif else begin
        fit = mpfitpeak(ygap,fgap,param,/moffat,nterms=6)
     endelse
;    For points that fall below 0, set them to 0
     ineg = where(fit lt 0,ctneg)
     if ctneg gt 0 then fit[ineg] = 0

     if doplot then begin
        cgplot,ygap,fgap,xran=[0,sx],yran=[min(fgap),max(fgap)],/xsty,/ysty
        cgplot,ygap,fit,/overplot,color='Blue'
        cgplot,ysec,xsec,/overplot,color='Red'
        wait,2
     endif

     yscat = moffat(ysec,param)
;    For points that fall below 0, set them to 0
     ineg = where(yscat lt 0,ctneg)
     if ctneg gt 0 then yscat[ineg] = 0
     yscatbin = transpose(rebin(yscat,sy,i2-i1+1))
     sci[i1:i2,*] -= yscatbin

  endfor

; Write output file
  mwrfits,phu,outfile,hdr_phu,/create,/silent
  mwrfits,sci,outfile,hdr_sci,/silent
  mwrfits,var,outfile,hdr_var,/silent
  mwrfits,dq,outfile,hdr_dq,/silent
  mwrfits,mdf,outfile,hdr_mdf,/silent

  sci = 0
  var = 0
  dq = 0

end
