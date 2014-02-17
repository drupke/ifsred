; docformat = 'rst'
;
;+
;
; Fit object centroids as a function of wavelength using IFSR_PEAK. Compute 
; differential atmospheric refraction correction based on object centroids. 
; Apply DAR correction to data cube.
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Postscript plots of centroid coordinates vs. wavelength and DAR corrected 
;    data cube. First column and row in plot are (1,1). Output plots have same 
;    filename as input cube, with "_dar.eps" at the end. Output cube has same 
;    filename as input cube, with "d" prepended, and has new header parameters
;    XPEAK and YPEAK, which are the peak determined at the reference wavelength.
;
; :Params:
;    incube: in, required, type=string
;      Path and filename for input data cube. Cube should contain data, 
;      variance, and DQ planes in extensions 1-3.
;    xord: in, required, type=integer
;      Order of polynomial fit to centroid columns vs. wavelength.
;    yord: in, required, type=integer
;      Order of polynomial fit to centroid rows vs. wavelength.
;
; :Keywords:
;    circ: in, optional, type=byte
;      Enforce fitting of a Gaussian symmetric in x and y.
;    dlam: in, optional, type=double, default=20
;      Number of Angstroms to sum for each wavelength bin in determining DAR.
;    dqthresh: in, optional, type=double, default=1
;      Threshold in DQ plane for assigning a pixel good/bad status. Greater 
;      than dqthresh means bad, less than or equal to dqthresh is good. Bad
;      pixels are set to DQ = 1.
;    lamcuts: in, optional, type=dblarr(2\,nregions)
;      Lower and upper limits of wavelength regions to ignore in determining 
;      DAR. E.g., [[5000,5005],[6800,6850],...]
;    lamref: in, optional, type=double, default=6300
;      Wavelength for which correction is set to 0.
;    nodarcor: in, optional, type=byte
;      Do not apply DAR correction to data cube.
;    sigcut: in, optional, type=double, default=3
;      Threshold in sigma for rejecting points in fit to centroids vs. 
;      wavelength.
;    xranfit: in, optional, type=dblarr(2)
;      Range for fitting centroid, in cube columns. Default is full range. First
;      column is 1.
;    xranplot: in, optional, type=dblarr(2)
;      Range for plotting centroid, in cube columns. Default is full range.
;      First column is 1.
;    yranfit: in, optional, type=dblarr(2)
;      Range for fitting centroid, in cube rows. Default is full range. First
;      row is 1.
;    yranplot: in, optional, type=dblarr(2)
;      Range for plotting centroid, in cube rows. Default is full range. First
;      row is 1.
;
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
;      2014feb04, DSNR, complete rewrite for ease of use/customization;
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
pro ifsr_dar,incube,xord,yord,circ=circ,dlam=dlam,$
             dqthresh=dqthresh,$
             lamcuts=lamcuts,lamref=lamref,$
             nodarcor=nodarcor,sigcut=sigcut,$
             xranplot=xranplot,xranfit=xranfit,$
             yranplot=yranplot,yranfit=yranfit


  if ~ keyword_set(dqthresh) then dqthresh = 1d
  if ~ keyword_set(dlam) then dlam = 20d
  if ~ keyword_set(lamref) then lamref=6300d
  if ~ keyword_set(circ) then circular=0 else circular=1
  ; Threshold in sigma for rejecting points
  if ~ keyword_set(sigcut) then sigcut = 3

; Load data cube
  header=1
  cube = ifsf_readcube(incube,header=header,/quiet)

; Set all DQs to 1 or 0
  gdind = where(cube.dq lt dqthresh,ctgdind)
  bdind = where(cube.dq ge dqthresh,ctbdind)
  cube.dq[gdind] = 0
  cube.dq[bdind] = 1

; Flux summation steps
  dz = round(dlam/cube.cdelt)
  npts = floor(double(cube.nz)/double(dz))

; Initialize outputs
  fc = dblarr(npts)
  xc = dblarr(npts)
  yc = dblarr(npts)
  wc = dblarr(npts)
  flag = intarr(npts)

;
;; Compute centroids
;
  lc_size = size(lamcuts)
; Make sure lamcuts array is 2-dimensional
  if lc_size[0] eq 1 then begin
     lamcuts = reform(lamcuts,lc_size[1],1)
     lc_size = size(lamcuts)
  endif
  iold=0
  j=0
; Cycle through wavelength steps
  for i=dz,cube.nz-1,dz do begin

;    Reject flagged wavelength regions
     for k=0,lc_size[2]-1 do begin
        if ((lamcuts[0,k] gt cube.wave[iold] AND $
             lamcuts[0,k] lt cube.wave[i]) OR $
            (lamcuts[1,k] gt cube.wave[iold] AND $
             lamcuts[1,k] lt cube.wave[i]) OR $
            (lamcuts[0,k] lt cube.wave[iold] AND $
             lamcuts[1,k] gt cube.wave[i])) then flag[j]=1
     endfor
     
;    Fit spatial peak
     peakpar = ifsr_peak(cube,[iold,i],circ=circular,/indrange,/quiet,$
                         xranfit=xranfit,yranfit=yranfit)
     xc[j]=peakpar[0]
     yc[j]=peakpar[1]
     wc[j]=(cube.wave[iold]+cube.wave[i])/2d
     if xc[j] lt 0 OR yc[j] lt 0 then flag[j]=1


     iold=i
     j++

  endfor

;
; Fit polynomial in wavelength space to centroid arrays
;


; X-direction

; Reject flagged points
  igx = where(~ flag,ngx)
  ibx = where(flag)
  wcgx=wc[igx]
  xcg=xc[igx]
  xcoeff=mpfitfun('poly',wcgx,xcg,0,dblarr(xord),$
                  weights=dblarr(ngx)+1d,yfit=xcfit,/quiet)
  if xord eq 0 then xcfit=dblarr(ngx)+xcfit

; Reject points by sigma*RMS deviation from fit
  xcg_rms=sqrt(median((xcg-xcfit)^2d))
  xcg_dev=abs(xcg-xcfit)
  ig2x=where(xcg_dev le sigcut*xcg_rms,nggx)
  ib2x=where(xcg_dev gt sigcut*xcg_rms)
  wcggx = wcgx[ig2x]
  xcgg = xcg[ig2x]
  xcoeff=mpfitfun('poly',wcggx,xcgg,0,dblarr(xord),$
                  weights=dblarr(nggx)+1d,yfit=xcfit,/quiet)
  if xord eq 0 then xcfit=dblarr(nggx)+xcfit

; Y-direction

  igy = where(~ flag,ngy)
  iby = where(flag)
  wcgy=wc[igy]
  ycg=yc[igy]
  ycoeff=mpfitfun('poly',wcgy,ycg,0,dblarr(yord),$
                  weights=dblarr(ngy)+1d,yfit=ycfit,/quiet)
  if yord eq 0 then ycfit=dblarr(ngy)+ycfit
  ycg_rms=sqrt(median((ycg-ycfit)^2d))
  ycg_dev=abs(ycg-ycfit)
  ig2y=where(ycg_dev le sigcut*ycg_rms,nggy)
  ib2y=where(ycg_dev gt sigcut*ycg_rms)
  wcggy = wcgy[ig2y]
  ycgg = ycg[ig2y]
  ycoeff=mpfitfun('poly',wcggy,ycgg,0,dblarr(yord),$
                  weights=dblarr(nggy)+1d,yfit=ycfit,/quiet)
  if yord eq 0 then ycfit=dblarr(nggy)+ycfit

;
;; Plot centroids
;

; Plot ranges
  nzx = where(xc ne 0)
  nzy = where(yc ne 0)
  if ~ keyword_set(xranplot) then xranplot = [min(xc[nzx]),max(xc[nzx])]
  if ~ keyword_set(yranplot) then yranplot = [min(yc[nzy]),max(yc[nzy])]
  lran = [min(cube.wave),max(cube.wave)]

; Strip ".fits" from input data file
  fitspos = strpos(incube,'.fits')
  if fitspos ne -1 then incubeuse=strmid(incube,0,fitspos) $
  else incubeuse=incube
  
  cgps_open,incubeuse+'_dar.eps',charsize=1,default_thick=2,/encap,$
     /inches,xs=10,ys=5,/qui

  multiplot,[2,1],gap=0.05
  cgplot,[0],xran=lran,yran=xranplot,/xsty,/ysty,$
       xtit=textoidl('Wavelength (A)'),ytit='X (spaxels)'
  symbols,2,1.5
; Good points
  cgoplot,wcggx,xcgg,psym=8
; Rejected points
  if ibx[0] ne -1 then cgoplot,wc[ibx],xc[ibx],psym=8,color='Green'
  if ib2x[0] ne -1 then cgoplot,wcgx[ib2x],xcg[ib2x],psym=8,color='Green'
; Fit
  cgoplot,wcggx,xcfit,color='Red',thick=4
  multiplot,/doyaxis
  cgplot,[0],xran=lran,yran=yranplot,/xsty,/ysty,$
       xtit=textoidl('Wavelength (A)'),ytit='Y (spaxels)'
  cgoplot,wcggy,ycgg,psym=8
  if iby[0] ne -1 then cgoplot,wc[iby],yc[iby],psym=8,color='Green'
  if ib2y[0] ne -1 then cgoplot,wcgy[ib2y],ycg[ib2y],psym=8,color='Green'
  cgoplot,wcggy,ycfit,color='Red',thick=4

  multiplot,/reset
  cgps_close

  xcc = poly(lamref,xcoeff)
  ycc = poly(lamref,ycoeff)

;
;; Apply DAR
;
  if ~ keyword_set(nodarcor) then begin

     for i=0,cube.nz-1 do begin

        dx = poly(cube.wave[i],xcoeff) - xcc[0]
        dy = poly(cube.wave[i],ycoeff) - ycc[0]

        flux = cube.dat[*,*,i]
        fluxvar = cube.var[*,*,i]
        fluxdq = cube.dq[*,*,i]
        xarr = dindgen(cube.ncols)+dx[0]
        yarr = dindgen(cube.nrows)+dy[0]

        newflux = interpolate(flux,xarr,yarr,/grid)
        newfluxvar = interpolate(fluxvar,xarr,yarr,/grid)
        newfluxdq = interpolate(fluxdq,xarr,yarr,/grid)
        
        cube.dat[*,*,i] = newflux
        cube.var[*,*,i] = newfluxvar
        cube.dq[*,*,i] = newfluxdq

     endfor

;    Add peak locations to header. SXADDPAR increases the size of the header, 
;    which confuses the structure, so we have to define a new array to hold the 
;    new header.
     newheader_phu = header.phu
     sxaddpar,newheader_phu,'XPEAK',xcc[0],'Column of peak flux (pixels)'
     sxaddpar,newheader_phu,'YPEAK',ycc[0],'Row of peak flux (pixels)'

;    Correct for DQ interpolation
     zeroind = where(cube.dq le 0.01,ctzeroind)
     oneind = where(cube.dq gt 0.01,ctoneind)
     cube.dq[zeroind] = 0
     cube.dq[oneind] = 1

     ; Format output file
     incubesplit = strsplit(incube,'/',/extract,count=nsplit)
     if nsplit gt 1 then begin
       ;    case of pathing
       prefix=strjoin(incubesplit[0:nsplit-2],'/') + '/'
       ;    case of full path (leading slash)
       if strpos(incube,'/') eq 0 then prefix = '/'+prefix
     endif else prefix=''
     outfile=prefix+'d'+incubesplit[nsplit-1]

     writefits,outfile,cube.phu,newheader_phu
     writefits,outfile,cube.dat,header.dat,/append
     writefits,outfile,cube.var,header.var,/append
     writefits,outfile,cube.dq,header.dq,/append

  endif

end
