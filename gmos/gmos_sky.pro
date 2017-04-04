; docformat = 'rst'
;
;+
;
; Apply wavelength shift to individual spaxels as determined from a sky line.
; 
; 
; The 
;
; :Categories:
;    GMOS
;
; :Returns:
;    1. Maps of sky line shifts in velocity space, FWHM of sky lines, and flux 
;       of sky lines.
;    2. Sky subtracted data file. Routine prepends the letter "s" to the input 
;       filename to get the output file.
;
; :Params:
;    datfile: in, required, type=string
;      Data file, with path.
;    fitfile: in, required, type=string
;      Output file from IFSFA (.xdr).
;    skyline: in, required, type=string
;      Name of sky line used for fit (as specified in IFSF_LINELIST)
;
; :Keywords:
;    sigrej: in, optional, type=double, default=3.0
;      Sigma value used to reject outliers in velocity, velocity width, and 
;      flux. Outliers are replaced with the mean value.
;    noshift: in, optional, type=byte
;      Do not apply velocity shifts.
;    plotlab: in, optional, type=byte
;      Add index labels to spaxels in plots
;    skyaps: in, optional, type=dblarr(Naps)
;      Option to choose sky apertures manually, if not all sky apertures are
;      to used in subtraction.
;    dqval: in, optional, type=double, default=1
;      Value for output DQ. NOTE: GFCUBE assumes 8!
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
;      2011jun02, DSNR, fixed what may be a bug ... seemed to be ignoring 
;                       spectra with negative velocity shifts.
;      2014jan31, DSNR, complete rewrite for ease of use/customization; 
;                       added detailed documentation;
;                       automated stitching of sky and data in plots 
;      2015jan05, DSNR, added option for spaxel labels in plots; revamped
;                       output stats
;      2015aug12, DSNR, added option to specify sky apertures to use in subtraction
;      2016jan28, DSNR, mods for data taken thru blue slit only
;      2017jan10, DSNR, DQVAL par
;
; :Copyright:
;    Copyright (C) 2014-2016 David S. N. Rupke
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
pro gmos_sky,datfile,fitfile,skyline,sigrej=sigrej,noshift=noshift,$
             plotlab=plotlab,skyaps=skyaps,dqval=dqval

  if keyword_set(plotlab) then dolab=1b else dolab=0b
  if ~ keyword_set(dqval) then dqval=1d

  fwhmtosig = 2d*sqrt(2d*alog(2d))

; Get data and fit results
  ftab_ext,datfile,[1,2,3,5],ap_mdf,x,y,beam
  restore,file=fitfile

; Strip ".fits" from input data file
  fitspos = strpos(datfile,'.fits')
  if fitspos ne -1 then datfileuse=strmid(datfile,0,fitspos) $
  else datfileuse=datfile

; Open file for stats
  openw,lunstats,datfileuse+'_skyline_stats.txt',/get_lun

; Get wavelength of sky line used in fit
  linelist = ifsf_linelist(skyline)
  
; Determine number of apertures and number of GMOS "slits" used (1 or 2)
;  naps = max(ap_mdf)
  naps = max(ap_mdf)-min(ap_mdf)+1
  if naps gt 750 then nslits=2 else nslits=1
  
; Populate arrays
  goodap = where(beam ne -1,ngoodaps)
  flux_all = dblarr(naps)
  vel_all = dblarr(naps)
  fwhm_all = dblarr(naps)
  flux_all[goodap] = emlflx['fc1',skyline,*,0]
  vel_all[goodap] = (emlwav['c1',skyline,*,0]/linelist[skyline]-1d)*299792d
  fwhm_all[goodap] = emlsig['c1',skyline,*,0]*fwhmtosig

; Subsets

; 1 slit
  if nslits eq 1 then begin
;   red slit
    if min(ap_mdf) lt 750 then i1 = where(ap_mdf le 750) $
;   blue slit
    else i1 = where(ap_mdf gt 750)
  endif else begin
; 2 slits
    i2 = where(ap_mdf gt 750)
  endelse

; spaxels with non-zero sky line fluxes
  nonzero = where(flux_all ne 0)
  i1nz = cmset_op(nonzero,'AND',i1)
  if nslits eq 2 then i2nz = cmset_op(nonzero,'AND',i2)

; Join sky and data sections in the x direction.
  xmean = mean(x)
  ihi = where(x gt xmean,cthi)
  ilo = where(x lt xmean,ctlo)
  ; which section is the sky?
  if cthi gt ctlo then begin
    isky = ilo
    isci = ihi
  endif else begin
    isky = ihi
    isci = ilo
  endelse
  if keyword_set(skyaps) then isky = skyaps
  xjoinlo = max(x[ilo])
  xdiff = xjoinlo - max(x[where(x lt xmean AND x ne xjoinlo)])
  xjoinhi = min(x[ihi])
  x[ilo] += xjoinhi - xjoinlo - 2d*xdiff

; Reject outliers
  if ~keyword_set(sigrej) then sigrej = 3d
  medvel1 = median(vel_all[i1nz])
  medfwhm1 = median(fwhm_all[i1nz])
  medflux1 = median(flux_all[i1nz])
  rmsvel1 = sqrt(median((vel_all[i1nz]-medvel1)^2d))
  rmsfwhm1 = sqrt(median((fwhm_all[i1nz]-medfwhm1)^2d))
  rmsflux1 = sqrt(median((flux_all[i1nz]-medflux1)^2d))
  ibad1 = cmset_op(i1nz,'AND',where(abs(vel_all-medvel1) $
                                    gt sigrej*rmsvel1),count=nbad1)
  vel_all[ibad1] = medvel1
  fwhm_all[ibad1] = medfwhm1
  flux_all[ibad1] = medflux1
  printf,lunstats,'------'
  printf,lunstats,'SLIT 1'
  printf,lunstats,'------'
  printf,lunstats,'VEL  median / RMS / med/RMS: ',$
                  medvel1,rmsvel1,rmsvel1/medvel1,$
                  format='(A-30,I10,D10.1,D10.3)'
  printf,lunstats,'FWHM median / RMS / med/RMS: ',$
                  medfwhm1,rmsfwhm1,rmsfwhm1/medfwhm1,$
                  format='(A-30,I10,D10.1,D10.3)'
  printf,lunstats,'FLUX median / RMS / med/RMS: ',$
                  medflux1,rmsflux1,rmsflux1/medflux1,$
                  format='(A-30,E10.2,E10.2,D10.3)'
  printf,lunstats,'No. pts. replaced: ',nbad1,format='(A-25,I0)'

  if nslits eq 2 then begin
     medvel2 = median(vel_all[i2nz])
     medfwhm2 = median(fwhm_all[i2nz])
     medflux2 = median(flux_all[i2nz])
     rmsvel2 = sqrt(median((vel_all[i2nz]-medvel2)^2d))
     rmsfwhm2 = sqrt(median((fwhm_all[i2nz]-medfwhm2)^2d))
     rmsflux2 = sqrt(median((flux_all[i2nz]-medflux2)^2d))
     ibad2 = cmset_op(i2nz,'AND',where(abs(vel_all-medvel2) $
                                       gt sigrej*rmsvel2),count=nbad2)
     vel_all[ibad2] = medvel2
     fwhm_all[ibad2] = medfwhm2
     flux_all[ibad2] = medflux2
     printf,lunstats,'------'
     printf,lunstats,'SLIT 2'
     printf,lunstats,'------'
     printf,lunstats,'VEL  median / RMS / med/RMS: ',$
                     medvel2,rmsvel2,rmsvel2/medvel2,$
                     format='(A-30,I10,D10.1,D10.3)'
     printf,lunstats,'FWHM median / RMS / med/RMS: ',$
                     medfwhm2,rmsfwhm2,rmsfwhm2/medfwhm2,$
                     format='(A-30,I10,D10.1,D10.3)'
     printf,lunstats,'FLUX median / RMS / med/RMS: ',$
                     medflux2,rmsflux2,rmsflux2/medflux2,$
                     format='(A-30,E10.2,E10.2,D10.3)'
  endif

; Sky vs. Sci stats
  printf,lunstats,'-----------'
  printf,lunstats,'SKY VS. SCI'
  printf,lunstats,'-----------'
  inzsci = cmset_op(nonzero,'AND',isci)
  medvelsci = median(vel_all[inzsci])
  medfwhmsci = median(fwhm_all[inzsci])
  medfluxsci = median(flux_all[inzsci])
  inzsky = cmset_op(nonzero,'AND',isky)
  medvelsky = median(vel_all[inzsky])
  medfwhmsky = median(fwhm_all[inzsky])
  medfluxsky = median(flux_all[inzsky])
  printf,lunstats,'VEL  sky / sci / frac diff: ',$
                  medvelsky,medvelsci,$
                  abs(medvelsky-medvelsci)/((medvelsky+medvelsci)/2d),$
                  format='(A-30,D10.1,D10.1,D10.3)'
  printf,lunstats,'FWHM  sky / sci / frac diff: ',$
                  medfwhmsky,medfwhmsci,$
                  abs(medfwhmsky-medfwhmsci)/((medfwhmsky+medfwhmsci)/2d),$
                  format='(A-30,D10.1,D10.1,D10.3)'
  printf,lunstats,'FLUX  sky / sci / frac diff: ',$
                  medfluxsky,medfluxsci,$
                  abs(medfluxsky-medfluxsci)/((medfluxsky+medfluxsci)/2d),$
                  format='(A-30,E10.2,E10.2,D10.3)'


; Set ranges and plot
  zran = [min(fwhm_all[nonzero]),max(fwhm_all[nonzero])]
  gmos_maps_fiber,nslits,fwhm_all[goodap],x[goodap],y[goodap],$
                  datfileuse+'_skyline_fwhm',zran=zran,dolab=dolab
  zran = [min(vel_all[nonzero]),max(vel_all[nonzero])]
  gmos_maps_fiber,nslits,vel_all[goodap],x[goodap],y[goodap],$
                  datfileuse+'_skyline_vel',cbform='(I0)',$
                  zran=zran,dolab=dolab
  zran = [min(flux_all[nonzero]),max(flux_all[nonzero])]
  gmos_maps_fiber,nslits,flux_all[goodap],x[goodap],y[goodap],$
                  datfileuse+'_skyline_flux',cbform='(E0.1)',$
                  zran=zran,dolab=dolab

; Correct line velocities

; Read data
  header=1
  cube=ifsf_readcube(datfile,header=header,/oned,dat=2,var=3,dq=4,/quiet)

  if not keyword_set(noshift) then begin
     for i=0,ngoodaps-1 do begin
        if vel_all[i] ne 0 then begin
           xarr = dindgen(cube.nz) + $
                  (vel_all[i]/299792d)*linelist[skyline]/cube.crpix
           newdat = interpolate(cube.dat[*,i],xarr)
           newvar = interpolate(cube.var[*,i],xarr)
           newdq = interpolate(cube.dq[*,i],xarr)
           zeroind = where(newdq le 0.01,ctzeroind)
           oneind = where(newdq gt 0.01,ctoneind)
           newdq[zeroind] = 0
           newdq[oneind] = dqval
           cube.dat[*,i] = newdat
           cube.var[*,i] = newvar
           cube.dq[*,i] = newdq
        endif
     endfor
  endif

; Subtract sky
  skyspec = median(cube.dat[*,isky],dim=2,/double)
  skymat = rebin(skyspec,cube.nz,cube.ncols)
  cube.dat -= skymat

; Format output file
  datfilesplit = strsplit(datfile,'/',/extract,count=nsplit)
  if nsplit gt 1 then begin
;    case of pathing
     prefix=strjoin(datfilesplit[0:nsplit-2],'/') + '/'
;    case of full path (leading slash)
     if strpos(datfile,'/') eq 0 then prefix = '/'+prefix
  endif else prefix=''
  outfile=prefix+'s'+datfilesplit[nsplit-1]

; Write output
  writefits,outfile,cube.phu,header.phu
  mdf = readfits(datfile,header_mdf,ext=1,/silent)
  writefits,outfile,mdf,header_mdf,/append
  writefits,outfile,cube.dat,header.dat,/append
  writefits,outfile,cube.var,header.var,/append
  writefits,outfile,cube.dq,header.dq,/append
  
  free_lun,lunstats

end
