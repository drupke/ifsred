; docformat = 'rst'
;
;+ 
; 
; Takes exposures ready for processing with NIFCUBE and co-aligns them using a 
; peak-finding routine, computes and applies a DAR correction, and interpolates
; the data and variance into data cubes.
;
; :Categories:
;    IFSRED/NIFS
;
; :Returns:
;    NIFS data cubes.
;
; :Params:
;    inlist: in, required, type=string
;      List of input FITS files.
;    outlist: in, required, type=string
;      List of output FITS files.
;    outpsf: in, required, type=string
;      File of peak fitting results.
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
;      2013feb26, DSNR, copied from GMOS routine
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
pro nifs_align,inlist,outlist,outpsf

;--------------------------------------------------------------------
; User inputs
;--------------------------------------------------------------------

; input size of spaxels in arcseconds
  dx_as_in = 0.1035d
  dy_as_in = 0.0445d
; output size of spaxels in arcseconds
  dx_as_out = dx_as_in
  dy_as_out = dx_as_in
; reference spectral coordinates
  crvalref = 18970.7d
  cdeltref = 2.135355d
; approximate RA and dec (from NED) in decimal degrees
  ra = 56.8736769d
  dec = 194.0593079d
; separation between points for DAR correction
  dwave = 200d
; geometry of input and output
  ncols_in = 29
  nrows_in = 69
  ncols_out = 45
  nrows_out = 45
  nz_in = 2040
  nz_out = 2040
; orders of polynomial fits to DAR (position vs. wavelength)
  xord = 2
  yord = 2

; Output wavelength solution     
  crpixref = 1
  pix = dindgen(nz_in)
  waveref = crvalref + cdeltref*(pix+1-crpixref) 

;--------------------------------------------------------------------
; Load data
;--------------------------------------------------------------------

; input arrays
  cube_dat_in = dblarr(ncols_in,nrows_in,nz_in)
  cube_var_in = dblarr(ncols_in,nrows_in,nz_in)
  cube_dq_in = dblarr(ncols_in,nrows_in,nz_in)

; Cycle through files
  openr,lun,inlist,/get_lun
  openr,lunout,outlist,/get_lun
  openw,lunpsf,outpsf,/get_lun
  ifile = 0
  infile=''
  outfile=''
  while (~ EOF(lun)) do begin

     readf,lun,infile
     readf,lunout,outfile

     dumy = readfits(infile,phu,ext=0,/silent)
     for i = 0,ncols_in-1 do begin
        cube_dat_in[i,*,*] = $
           transpose(readfits(infile,hdr,ext=3*i+1,/silent))
        cube_var_in[i,*,*] = $
           transpose(readfits(infile,hdr_var,ext=3*i+2,/silent))
        cube_dq_in[i,*,*] = $
           transpose(readfits(infile,hdr_dq,ext=3*i+3,/silent))
     endfor
;     cube_var_in = sqrt(abs(cube_var_in))

; Flip around y axis
     cube_dat_in = reverse(cube_dat_in,1)
     cube_var_in = reverse(cube_var_in,1)
     cube_dq_in = reverse(cube_dq_in,1)

;--------------------------------------------------------------------
; Find centroids
;--------------------------------------------------------------------

; Fit ranges
; Collapse image to estimate center
     img = total(cube_dat_in[*,*,99:nz_in-100],3)
     imgvar = total(cube_var_in[*,*,99:nz_in-100],3)
     pkval = max(img,ipkval)
     xypkval = array_indices(img,ipkval)
     xranfit = [xypkval[0]-5,xypkval[0]+5]
     yranfit = [xypkval[1]-10,xypkval[1]+10]

; Now get a better handle on the center ...

     xpkfit = dindgen(xranfit[1]-xranfit[0]+1)+xranfit[0]
     ypkfit = dindgen(yranfit[1]-yranfit[0]+1)+yranfit[0]
     dat_xpkfit = total(img[xranfit[0]:xranfit[1],yranfit[0]:yranfit[1]],2)
     dat_ypkfit = total(img[xranfit[0]:xranfit[1],yranfit[0]:yranfit[1]],1)
     var_xpkfit = total(imgvar[xranfit[0]:xranfit[1],yranfit[0]:yranfit[1]],2)
     var_ypkfit = total(imgvar[xranfit[0]:xranfit[1],yranfit[0]:yranfit[1]],1)

     fit_xpkfit = mpfitpeak(xpkfit,dat_xpkfit,xparam,/lorentz,$
                            weights=1d/var_xpkfit,nterms=4)
     fit_ypkfit = mpfitpeak(ypkfit,dat_ypkfit,yparam,/lorentz,$
                            weights=1d/var_ypkfit,nterms=4)
     pk_xpkfit = xparam[0]
     pk_ypkfit = yparam[0]
     ;; fw_xpkfit = xparam[2] * 2.35 * dx_as_in
     ;; fw_ypkfit = yparam[2] * 2.35 * dy_as_in
     fw_xpkfit = xparam[2] * 2 * dx_as_in
     fw_ypkfit = yparam[2] * 2 * dy_as_in

;    Print peak fits to the screen and to a file
     print,ifile+1,infile,pkval,fw_xpkfit,fw_ypkfit,$
           mean([fw_xpkfit,fw_ypkfit]),$
           format='(I-3,A0,I10,D6.3,D6.3,D6.3)'
     printf,lunpsf,ifile+1,infile,pkval,fw_xpkfit,fw_ypkfit,$
           mean([fw_xpkfit,fw_ypkfit]),$
           format='(I-3,A0,I10,D6.3,D6.3,D6.3)'



;--------------------------------------------------------------------
; Plot centroid fit
;--------------------------------------------------------------------

     cleanplot,/silent
     set_plot,'ps'
     device,/inches,xsize=5,ysize=5,xoffset=0.5,yoffset=0.5
     device,/encapsulated,/color
     loadct,0,/silent
     !P.thick=2
     !P.charsize=1
     !P.charthick=2

     flran_xpkfit = [min([dat_xpkfit,fit_xpkfit]),$
                     max([dat_xpkfit,fit_xpkfit])]
     flran_ypkfit = [min([dat_ypkfit,fit_ypkfit]),$
                     max([dat_ypkfit,fit_ypkfit])]

     tmpstr = strsplit(infile,'.',/extract)
     device,filename=tmpstr[0]+'_pk.eps'

     cgplot,[0],xran=xranfit,yran=flran_xpkfit,/xsty,/ysty,$
            xtit='X (spaxels)',ytit='flux',$
            layout=[1,2,1]
     cgplot,xpkfit,dat_xpkfit,/over,thick=4
     cgplot,xpkfit,fit_xpkfit,/over,color='Red',thick=4
     loadct,13,/silent
     cgplot,[0],xran=yranfit,yran=flran_ypkfit,/xsty,/ysty,$
            xtit='Y (spaxels)',ytit='flux',$
            layout=[1,2,2]
     cgplot,ypkfit,dat_ypkfit,/over,thick=4
     cgplot,ypkfit,fit_ypkfit,/over,color='Red',thick=4
     loadct,13,/silent

     device,/close_file
     multiplot,/reset

; Compute wavelengths
     crval = double(sxpar(hdr,'CRVAL1',/silent))
     crpix = double(sxpar(hdr,'CRPIX1',/silent))
     cdelt = double(sxpar(hdr,'CDELT1',/silent))
     wave = crval + cdelt*(pix+1-crpix) 

; Flux summation steps
     dz = round(dwave/cdelt)
     npts = floor(double(nz_in)/double(dz))

; Initialize outputs
     xc = dblarr(npts)
     yc = dblarr(npts)
     wc = dblarr(npts)
     flag = intarr(npts)
; go ahead and flag first and last points
     flag[0] = 1
     flag[npts-1] = 1

; Compute centroids
     iold=0
     j=0
; Cycle through wavelength steps
     for i=dz,nz_in-1,dz do begin
;    Subsets of data to find peak in
        datsub = cube_dat_in[xranfit[0]:xranfit[1],yranfit[0]:yranfit[1],iold:i]
        varsub = cube_var_in[xranfit[0]:xranfit[1],yranfit[0]:yranfit[1],iold:i]
        dqsub = cube_dq_in[xranfit[0]:xranfit[1],yranfit[0]:yranfit[1],iold:i]
;    Zero out DQed data
        ;; ibad = where(dqsub eq 1)
        ;; datsub[ibad] = 0
        ;; varsub[ibad] = 0
;    Sum data in wavelength direction
        flux = total(datsub,3,/double)
        fluxvar = total(varsub,3,/double)
;    Normalize by number of wavelength pixels summed
        fluxct = (i-iold+1)  ;- total(dqsub,3,/double)
        flux /= fluxct
        fluxvar /= fluxct
;    Flag spaxels with no data
        inoflux = where(fluxct eq 0,ctnoflux)
        if ctnoflux gt 0 then begin
           flux[inoflux]=0d
           fluxvar[inoflux]=1d99
        endif
;    Fit spatial peak. I think this returns 0-index centers.
        yfit = mpfit2dpeak(flux,param,perror=perror,/lorentz,$
                           weights=1d/fluxvar)
        xc[j]=param[4]+xranfit[0]
        yc[j]=param[5]+yranfit[0]
        wc[j]=(wave[iold]+wave[i])/2d
        if xc[j] lt 0 OR yc[j] lt 0 then flag[j]=1
        iold=i
        j++
     endfor

;--------------------------------------------------------------------
; Fit polynomial in wavelength space to DAR correction
;--------------------------------------------------------------------

; Threshold in sigma for rejecting points
     badthresh = 5

; X-direction
; Reject flagged points
     igx = where(~ flag,ngx)
     ibx = where(flag)
     wcgx=wc[igx]
     xcg=xc[igx]
; First fit
     xcoeff=mpfitfun('poly',wcgx,xcg,0,dblarr(xord),$
                     weights=dblarr(ngx)+1d,yfit=xcfit,/quiet)
; Reject points by sigma*RMS deviation from fit
     xcg_rms=sqrt(median((xcg-xcfit)^2d))
     xcg_dev=abs(xcg-xcfit)
     ig2x=where(xcg_dev le badthresh*xcg_rms,nggx)
     ib2x=where(xcg_dev gt badthresh*xcg_rms)
     wcggx = wcgx[ig2x]
     xcgg = xcg[ig2x]
; Second fit, minus rejected points
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
     ig2y=where(ycg_dev le badthresh*ycg_rms,nggy)
     ib2y=where(ycg_dev gt badthresh*ycg_rms)
     wcggy = wcgy[ig2y]
     ycgg = ycg[ig2y]
     ycoeff=mpfitfun('poly',wcggy,ycgg,0,dblarr(yord),$
                     weights=dblarr(nggy)+1d,yfit=ycfit,/quiet)
     if yord eq 0 then ycfit=dblarr(nggy)+ycfit

;--------------------------------------------------------------------
; Plot DAR fit
;--------------------------------------------------------------------

     cleanplot,/silent
     set_plot,'ps'
     device,/inches,xsize=10,ysize=5,xoffset=0.5,yoffset=0.5
     device,/encapsulated,/color
     loadct,0,/silent
     !P.thick=2
     !P.charsize=1
     !P.charthick=2

     nzx = where(xc ne 0 and flag ne 1)
     nzy = where(yc ne 0 and flag ne 1)
     xmin = min(xc[nzx])
     xmax = max(xc[nzx])
     dxmm = xmax - xmin
     ymin = min(yc[nzy])
     ymax = max(yc[nzy])
     dymm = ymax - ymin
     xranplot = [xmin - 0.5*dxmm,xmax + 0.5*dxmm]
     yranplot = [ymin - 0.5*dymm,ymax + 0.5*dymm]
     lran = [min(wave),max(wave)]

     tmpstr = strsplit(infile,'.',/extract)
     device,filename=tmpstr[0]+'_dar.eps'

     multiplot,[2,1],gap=0.05
     cgplot,[0],xran=lran,yran=xranplot,/xsty,/ysty,$
            xtit='Wavelength ($\angstrom$)',ytit='X (spaxels)'
     symbols,2,1.5
; Good points
     oplot,wcggx,xcgg,psym=8
     loadct,13,/silent
; Rejected points
     if ibx[0] ne -1 then oplot,wc[ibx],xc[ibx],psym=8,color=125
     if ib2x[0] ne -1 then oplot,wcgx[ib2x],xcg[ib2x],psym=8,color=125
; Fit
     oplot,wcggx,xcfit,color=255,thick=4
     loadct,0,/silent
     multiplot,/doyaxis
     cgplot,[0],xran=lran,yran=yranplot,/xsty,/ysty,$
            xtit='Wavelength ($\angstrom$)',ytit='Y (spaxels)'
     oplot,wcggy,ycgg,psym=8
     loadct,13,/silent
     if iby[0] ne -1 then oplot,wc[iby],yc[iby],psym=8,color=125
     if ib2y[0] ne -1 then oplot,wcgy[ib2y],ycg[ib2y],psym=8,color=125
     oplot,wcggy,ycfit,color=255,thick=4
     loadct,0,/silent

     device,/close_file
     multiplot,/reset

;--------------------------------------------------------------------
; Interpolate spatially
;--------------------------------------------------------------------

;    x and y indices of center of input cube
     ixc_in = floor(ncols_in/2d)
     iyc_in = floor(nrows_in/2d)
;    x and y indices of center of output cube
     ixc_out = floor(ncols_out/2d)
     iyc_out = floor(nrows_out/2d)

     cube_dat_out = dindgen(nrows_out,ncols_out,nz_out)
     cube_var_out = dindgen(nrows_out,ncols_out,nz_out)
     cube_dq_out = dindgen(nrows_out,ncols_out,nz_out)

     for i=0,nz_in-1 do begin
;       x and y indices of center of slice
        ixc_slice = poly(wave[i],xcoeff)
        iyc_slice = poly(wave[i],ycoeff)
;       x and y arrays of locations of desired output points, in input
;       coordinates
        ;; xarr = (dindgen(ncols_out) - ixc_out)*dx_as_out/dx_as_in + ixc_slice
        ;; yarr = dindgen(nrows_out) - iyc_out + iyc_slice
        xarr = dindgen(ncols_out) - ixc_out + ixc_slice
        yarr = (dindgen(nrows_out) - iyc_out)*dy_as_out/dy_as_in + iyc_slice

        cube_dat_out[*,*,i] = $
           interpolate(cube_dat_in[*,*,i],xarr,yarr,/grid,missing=0d)
        cube_var_out[*,*,i] = $
           interpolate(cube_var_in[*,*,i],xarr,yarr,/grid,missing=0d)
        cube_dq_out[*,*,i] = $
           interpolate(cube_dq_in[*,*,i],xarr,yarr,/grid,missing=1d)

     endfor

;--------------------------------------------------------------------
; Interpolate spectrally
;--------------------------------------------------------------------

     zarr = dindgen(nz_out) + (waveref - wave)/cdeltref
     xarr = dindgen(ncols_out)
     yarr = dindgen(nrows_out)

     cube_dat_out = interpolate(cube_dat_out,xarr,yarr,zarr,/grid)
     cube_var_out = interpolate(cube_var_out,xarr,yarr,zarr,/grid)
     cube_dq_out = interpolate(cube_dq_out,xarr,yarr,zarr,/grid)

;--------------------------------------------------------------------
; Write back to disk
;--------------------------------------------------------------------

     sxaddpar,hdr,'EXTVER',1
     sxaddpar,hdr,'BITPIX',-64
     sxaddpar,hdr,'NAXIS',2
     sxaddpar,hdr,'PEAKVAL',pkval
     sxaddpar,hdr,'PEAKVALX',pk_xpkfit
     sxaddpar,hdr,'PEAKVALY',pk_ypkfit
     sxaddpar,hdr,'XFWHM',fw_xpkfit
     sxaddpar,hdr,'YFWHM',fw_ypkfit
     sxaddpar,hdr,'DISPAXIS',3
     sxaddpar,hdr,'CRPIX1',double(ixc_out+1)
     sxaddpar,hdr,'CRPIX2',double(iyc_out+1)
     sxaddpar,hdr,'CRPIX3',double(crpixref)
     sxaddpar,hdr,'CRVAL1',ra
     sxaddpar,hdr,'CRVAL2',dec
     sxaddpar,hdr,'CRVAL3',crvalref
     sxaddpar,hdr,'CD1_1',dx_as_out
     sxaddpar,hdr,'CD2_2',dy_as_out
     sxaddpar,hdr,'CD3_3',cdeltref
     sxaddpar,hdr,'LTM1_1',double(1)
     sxaddpar,hdr,'LTM2_2',double(1)
     sxaddpar,hdr,'LTM3_3',double(1)
     sxaddpar,hdr,'WAT0_001','system=world'
     sxaddpar,hdr,'WAT1_001','wtype=linear axtype=xi'
     sxaddpar,hdr,'WAT2_001','wtype=linear axtype=eta'
     sxaddpar,hdr,'WAT3_001','wtype=linear axtype=wave'
     sxdelpar,hdr,'WAXMAP01'
     sxdelpar,hdr,'CDELT1'
     sxdelpar,hdr,'CDELT2'

     hdr_var = hdr
     hdr_dq = hdr
     sxaddpar,hdr_var,'EXTNAME','VAR'
     sxaddpar,hdr_dq,'EXTNAME','DQ'
     
     zeroind = where(cube_dq_out le 0.1,ctzeroind)
     oneind = where(cube_dq_out gt 0.1,ctoneind)
     cube_dq_out[zeroind] = 0
     cube_dq_out[oneind] = 1
  
     writefits,outfile,dumy,phu
     writefits,outfile,cube_dat_out,hdr,/append
     writefits,outfile,cube_var_out,hdr_var,/append
     writefits,outfile,cube_dq_out,hdr_dq,/append
     
     ifile++
     
  endwhile

  free_lun,lun
  free_lun,lunout
  free_lun,lunpsf
  
  cube_dat_in = 0
  cube_dat_out = 0
  cube_var_in = 0
  cube_var_out = 0
  cube_dq_in = 0
  cube_dq_out = 0

end
