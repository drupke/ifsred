; docformat = 'rst'
;
;+
;
; Fit a 2D Gaussian to a data cube over a specified wavelength range using 
; MPFIT2DPEAK. The data is first summed over that range.
;
; :Categories:
;    IFSRED
;
; :Returns:
;   Array of peak column and row and half-widths. First column and row are 
;   (1,1).
;
; :Params:
;    cube: in, required, type=structure
;      Contains data cube and associated information, in the form output by 
;      IFSF_READCUBE.
;    lamrange: in, required, type=dblarr(2)
;      Lower and upper limits of wavelength regions to include in fit.
;    
; :Keywords:
;    circ: in, optional, type=byte
;      Enforce fitting of a Gaussian symmetric in x and y.
;    header: in, optional, type=structure
;      Header structure (output from IFSF_READCUBE) used for output file
;      if re-writing header with new XPEAK, YPEAK.
;    indrange: in, optional, type=byte
;      Give range for wavelength region to fit as indices rather
;      than actual wavelengths. Indices are put in lamrange parameter.
;    outfile: in, optional, type=string
;      Optional output filename if re-writing header with new XPEAK, YPEAK.
;    nophu: in, optional, type=byte
;      Added flag to indicate that the 0th extension contains the data, not 
;      the PHU.
;    quiet: in, optional, type=byte
;      Suppress output of centroids and half-widths to terminal.
;    xranfit: in, optional, type=dblarr(2)
;      Range for fitting centroid, in cube columns, in single-offset indices. 
;      Default is full range.
;    yranfit: in, optional, type=dblarr(2)
;      Range for fitting centroid, in cube rows, in single-offset indices. 
;      Default is full range.
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
;      2011may10, DSNR, improved rejection of outliers
;      2014feb06, DSNR, complete rewrite for ease of use/customization;
;                       added detailed documentation
;      2015may18, DSNR, added option to re-write header
;      2016sep12, DSNR, added NOPHU option
;      2020jan28, DSNR, added NaN rejection
;      2021dec07, DSNR, bugfix in treatment of extensions
;                       
; :Copyright:
;    Copyright (C) 2014--2021 David S. N. Rupke
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
function ifsr_peak,filename,lamrange,circ=circ,indrange=indrange,quiet=quiet,$
   xranfit=xranfit,yranfit=yranfit,noupdate=noupdate,$
   moffat=moffat,datext=datext,$
   varext=varext,dqext=dqext,funit=funit

   if ~keyword_set(datext) then datext=1
   if ~keyword_set(varext) then varext=2
   if ~keyword_set(dqext) then dqext=3
   if ~keyword_set(funit) then funit=1d
   if ~ keyword_set(circ) then circular=0 else circular=1

   header=1b
   cube = ifsf_readcube(filename,quiet=quiet,datext=datext,varext=varext,$
      dqext=dqext,header=header)

;  Indices of fit range in wavelength
   if ~ keyword_set(indrange) then begin
      indrange = intarr(2)
      indrange[0] = value_locate(cube.wave,lamrange[0])
      indrange[1] = value_locate(cube.wave,lamrange[1])
   endif else begin
      indrange = lamrange
   endelse
  
; Fit ranges in X and Y. Else part converts columns and rows that start at 1 
; to zero-offset indices.
   if ~keyword_set(xranfit) then xranfituse=[0,cube.ncols-1] $
   else xranfituse = xranfit - 1
   ncolsuse = xranfit[1]-xranfit[0]+1
   if ~keyword_set(yranfit) then yranfituse=[0,cube.nrows-1] $
   else yranfituse = yranfit - 1
   nrowsuse = yranfit[1]-yranfit[0]+1

   wavsub = cube.wave[indrange[0]:indrange[1]]
   datsub = cube.dat[xranfituse[0]:xranfituse[1],$
                     yranfituse[0]:yranfituse[1],$
                     indrange[0]:indrange[1]]
   varsub = cube.var[xranfituse[0]:xranfituse[1],$
                     yranfituse[0]:yranfituse[1],$
                     indrange[0]:indrange[1]]
   dqsub = cube.dq[xranfituse[0]:xranfituse[1],$
                   yranfituse[0]:yranfituse[1],$
                   indrange[0]:indrange[1]]
   nlamuse = n_elements(wavsub)

;  Zero out DQed or INF/NANed data
   ibad = where(dqsub eq 1, ctbad)
   if ctbad gt 0 then begin
      datsub[ibad] = 0d
      varsub[ibad] = 0d
   endif
   iinf = where(~ finite(datsub), ctinf)
   if ctinf gt 0 then begin
      datsub[iinf] = 0d
      varsub[iinf] = 0d
      dqsub[iinf] = 1
   endif

;  Sum over wavelength space
   dlam = wavsub[1:nlamuse-1]-wavsub[0:nlamuse-2]
   dlam = [dlam[0],dlam]
   dlamcube = rebin(reform(dlam,1,1,nlamuse),ncolsuse,nrowsuse,nlamuse)
   datsum = total(datsub*dlamcube,3,/double)*funit
   varsum = total(varsub*dlamcube^2d,3,/double)*funit^2d

;  Normalize to the same # of wavelength pixels summed at each spaxel
   fluxct = (indrange[1]-indrange[0]+1) - total(dqsub,3,/double)
   datsum /= fluxct/max(fluxct)
   varsum /= fluxct/max(fluxct)

;  Flag spaxels with no data
   inoflux = where(fluxct eq 0,ctnoflux)
   if ctnoflux gt 0 then begin
      datsum[inoflux]=0d
      varsum[inoflux]=1d99
   endif
  
   if keyword_set(moffat) then begin
      moffat=1b
;      npar=8
   endif else begin
      moffat=0b
;      npar=7
   endelse
;  parinfo = REPLICATE({fixed:0b, value:0d},npar)
;  if keyword_set(fixbaseline) then $
;    parinfo[0].fixed = 1b else parinfo[0].fixed=0b

   yfit = mpfit2dpeak(datsum,param,weights=1d/varsum,$
      circular=circular,moffat=moffat) ;,parinfo=parinfo)
;  Output parameters are central flux, sigmas or half widths, and centroids. +1 converts 
;  from 0-offset indices to cols and rows that start at 1.
   fc=param[1]
   xhw=param[2]
   yhw=param[3]
   xc=param[4]+xranfituse[0]+1
   yc=param[5]+yranfituse[0]+1
   if keyword_set(moffat) then begin
      xfwhm = 2d*xhw
      yfwhm = 2d*yhw
      gamma = param[7]
   endif else begin
      xfwhm = 2.35d*xhw
      yfwhm = 2.35d*yhw
   endelse

   if not keyword_set(quiet) then $
      print,'FPEAK/XPEAK/YPEAK/XFWHM/YFWHM = ',max(datsum),xc,yc,xfwhm,yfwhm,$
         format='(A0,E8.2,D8.2,D8.2,D6.2,D6.2)'

;
;; Re-write header with new XPEAK/YPEAK to new output file
;
   if not keyword_set(noupdate) then begin

     outfile=filename

;    Add peak locations to header. SXADDPAR increases the size of the header,
;    which confuses the structure, so we have to define a new array to hold the
;    new header.
     if datext eq 1 then newheader = header.phu $
     else newheader = header.dat
     sxaddpar,newheader,'FPEAK',fc,'Observed peak flux'
     sxaddpar,newheader,'XPEAK',xc,'Column of fitted peak (pixels)'
     sxaddpar,newheader,'YPEAK',yc,'Row of fitted peak (pixels)'
     sxaddpar,newheader,'XFWHM',xfwhm,'Column FWHM of PSF fit (pixels)'
     sxaddpar,newheader,'YFWHM',yfwhm,'Row FWHM of PSF fit (pixels)'

     if datext eq 1 then begin
        writefits,outfile,cube.phu,newheader
        writefits,outfile,cube.dat,header.dat,/append
     endif else begin
        writefits,outfile,cube.dat,newheader
     endelse
     writefits,outfile,cube.var,header.var,/append
     writefits,outfile,cube.dq,header.dq,/append

  endif

  return,[xc,yc,xfwhm,yfwhm]


end
