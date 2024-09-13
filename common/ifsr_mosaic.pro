; docformat = 'rst'
;
;+
;
; Mosaic data cubes for individual science exposures into a single
; data cube.  Accounts for shifts due to dithering, etc., using peaks
; determined with IFSR_DAR and written to the image headers as XPEAK
; and YPEAK. The data, variance, and data quality are interpolated to
; a common grid using INTERPOLATE and bilinear interpolation, and
; median-combined. The DQ plane is applied after interpolation but
; before median-combining.
;
; The new grid has the common reference coordinate (XPEAK, YPEAK) centered on
; a spaxel. The 
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Combined data cube. Header of cube is the same as the header of
;    the first exposure in the file list.
;    
;    Note that data is median-combined, so that output flux is not summed. 
;    Thus expected input flux is per unit time.
;    
; :Params:
;    infiles: in, required, type=strarr
;      Filenames for data cubes of individual exposures.
;    outfile: in, required, type=string
;      Path and filename of output file.
;
; :Keywords:
;    crpix: in, optional, type=dblarr(2)
;    crval: in, optional, type=dblarr(2)
;      Directly set spatial reference point of WCS for output image.
;    indir: in, optional, type=string
;      Directory where input files are located.
;    nocenter: in, optional, type=byte
;      Do not center the peak in the output data cube. Default is to
;      center it. Note that the this option presently may not result in the 
;      full dimensions of the data being output, for reasons TBD.
;    nophu: in, optional, type=byte
;      Added flag to indicate that the 0th extension contains the data, not
;      the PHU.
;    osiris: in, optional, type=byte
;    refexp: in, optional, type=integer
;      Reference exposure for header and spatial WCS
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
;      2010jun02, DSNR, created
;      2014feb06, DSNR, complete rewrite for ease of use/customization;
;                       added detailed documentation
;      2014aprXY, DSNR, fixed bug in call to SXPAR
;      2015aug03, DSNR, note made about NOCENTER keyword; needs to be fixed
;      2016sep12, DSNR, added NOPHU option
;      2018sep28, DSNR, small bug fixes for doing large mosaics where
;                       source not centered
;      2018nov04, DSNR, bug fix of variance computation. Was overestimating
;                       variance by factor of Nexp
;      2022oct25, DSNR, better commenting of algorithm. Added handling of 
;                       CRPIX, CRVAL, using reference exposure
;
; :Copyright:
;    Copyright (C) 2014--2022 David S. N. Rupke
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
pro ifsr_mosaic,infiles,outfile,indir=indir,nocenter=nocenter,nophu=nophu,$
   refexp=refexp,crpix=crpix,crval=crval,quiet=quiet

; Value for data outside the bounds of interpolation
  noval=-999d

  if not keyword_set(indir) then indir=''
  if not keyword_set(refexp) then refexp=1
  if not keyword_set(quiet) then quiet=1b

  nexp = n_elements(infiles)
  xcc = dblarr(nexp)
  ycc = dblarr(nexp)
  xccd = dblarr(nexp)
  yccd = dblarr(nexp)

; Loading data

  print,'IFSR_MOSAIC: Loading data.'

  for i=0,nexp-1 do begin

     header=1
     if keyword_set(nophu) then begin
        datext=-1
        varext=1
        dqext=2
     endif else begin
        datext=1
        varext=2
        dqext=3
     endelse
     cube = ifsf_readcube(indir+infiles[i],header=header,quiet=quiet,$
                          datext=datext,varext=varext,dqext=dqext)

;    Set output header to reference exposure
     if refexp eq i+1 then outheader=header

     if ~ keyword_set(nophu) then begin
        pkheader = header.phu
        appenddat=1
     endif else begin
        pkheader = header.dat
        appenddat=0
     endelse
     
;    Fitted galaxy center in zero-offset coordinates
     xcc[i] = sxpar(pkheader,'XPEAK',silent=quiet,count=ctx) - 1
     ycc[i] = sxpar(pkheader,'YPEAK',silent=quiet,count=cty) - 1
     if ~ keyword_set(nophu) AND (ctx eq 0 OR cty eq 0) then begin
;       Check science header if it's missing from PHU
        xcc[i] = sxpar(header.dat,'XPEAK','exposure'+string(i+1,format='(I0)'),$
                       /silent) - 1
        ycc[i] = sxpar(header.dat,'YPEAK','exposure'+string(i+1,format='(I0)'),$
                       /silent) - 1
     endif
;    Residual from next lowest integer pixel
     xccd[i] = xcc[i] - floor(xcc[i])
     yccd[i] = ycc[i] - floor(ycc[i])
     
;    Get spatial WCS information. Assume spaxel sizes all the same
     if refexp eq i+1 then begin
        crval1 = sxpar(pkheader,'CRVAL1',silent=quiet)
        crpix1 = sxpar(pkheader,'CRPIX1',silent=quiet)
        crval2 = sxpar(pkheader,'CRVAL2',silent=quiet)
        crpix2 = sxpar(pkheader,'CRPIX2',silent=quiet)
        crpix1 -= xcc[i]
        crpix2 -= ycc[i]
     endif

;    Initialize arrays to hold all exposures, as well as some other 
;    intermediate arrays
     if i eq 0 then begin
        nz0 = cube.nz
        crpix0 = cube.crpix
        crval0 = cube.crval
        cdelt0 = cube.cdelt
        datall = dblarr(cube.ncols,cube.nrows,cube.nz,nexp)
        varall = dblarr(cube.ncols,cube.nrows,cube.nz,nexp)
        dqall = dblarr(cube.ncols,cube.nrows,cube.nz,nexp)
        ; 2d arrays to hold x values and y values from original exposure,
        ; in zero-offset coordinates. These will be the coordinates of this
        ; exposure in the new grid system.
        xarr = rebin(dindgen(cube.ncols),cube.ncols,nexp)
        yarr = rebin(dindgen(cube.nrows),cube.nrows,nexp)
        ; Coordinates of the spaxels in this exposure with reference to peak
        xarrc = rebin(dindgen(cube.ncols),cube.ncols,nexp)
        yarrc = rebin(dindgen(cube.nrows),cube.nrows,nexp)
        ; This is the new grid coordinate system.
        xnewarr = rebin(indgen(cube.ncols),cube.ncols,nexp)
        ynewarr = rebin(indgen(cube.nrows),cube.nrows,nexp)
        ;
        newncols = cube.ncols
        newnrows = cube.nrows
     endif

;    Put everything on a common wavelength solution, if it is not already.
;    This bit assumes that the differences between wavelength solutions are found
;    in the initial pixel CRPIX. This may not always apply.
     if cube.cdelt ne cdelt0 then print,'WARNING: Exposures have different '+$
        'dispersions. IFSR_MOSAIC does not account for this.'
     dcrpix = cube.crpix - crpix0
     dcrval = round((cube.crval - crval0)/cdelt0)
     dcrval_fracdiff = (cube.crval - crval0)/cdelt0 - dcrval
     ;ddcrval_fracdiff = abs(((cube.crval - crval0)/cdelt0 - dcrval)/dcrval)
     if dcrval_fracdiff ge 0.05d  then print,'NOTE: Exposure '+$
        string(i+1,format='(I0)')+' is offset in lambda zero-point by '+$
        'fractional difference of '+string(dcrval_fracdiff,format='(D0.2)')+$
        ' dispersion elements, w.r.t. the reference exposure. IFSR_MOSAIC'+$
        ' does not account for this.'
     dpix = dcrpix - dcrval
     nzuse = (cube.nz le nz0) ? cube.nz : nz0
     ncolsuse = (cube.ncols le newncols) ? cube.ncols : newncols
     nrowsuse = (cube.nrows le newnrows) ? cube.nrows : newnrows
     if dpix le 0 then begin
        datall[0:ncolsuse-1,0:nrowsuse-1,abs(dpix):nzuse-1,i] = $
           cube.dat[0:ncolsuse-1,0:nrowsuse-1,0:nzuse-1-abs(dpix)]
        varall[0:ncolsuse-1,0:nrowsuse-1,abs(dpix):nzuse-1,i] = $
           cube.var[0:ncolsuse-1,0:nrowsuse-1,0:nzuse-1-abs(dpix)]
        dqall[0:ncolsuse-1,0:nrowsuse-1,abs(dpix):nzuse-1,i] = $
           cube.dq[0:ncolsuse-1,0:nrowsuse-1,0:nzuse-1-abs(dpix)]
     endif else begin
        datall[0:ncolsuse-1,0:nrowsuse-1,0:nzuse-1-dpix,i] = $
           cube.dat[0:ncolsuse-1,0:nrowsuse-1,dpix:nzuse-1]
        varall[0:ncolsuse-1,0:nrowsuse-1,0:nzuse-1-dpix,i] = $
           cube.var[0:ncolsuse-1,0:nrowsuse-1,dpix:nzuse-1]
        dqall[0:ncolsuse-1,0:nrowsuse-1,0:nzuse-1-dpix,i] = $
           cube.dq[0:ncolsuse-1,0:nrowsuse-1,dpix:nzuse-1]
     endelse

     ; Add in non-integer part of galaxy center position to coordinates of this
     ; exposure in the new grid coordinates. In these coordinates, peak is
     ; at [XPEAK, YPEAK] - 1
     ; e.g., if [XPEAK, YPEAK] = [23.42,28.08], then in these coordinates
     ; it's [22.42, 27.08]
     ; But actual x values are:
     ; [0.42, ..., 21.42, 22.42, 23.42, ...]
     xarr[*,i] += xccd[i]
     yarr[*,i] += yccd[i]
     ; Subtract peak coordinate from coordinates of this spaxel --
     ; in these coordinates, peak is at [0,0].
     ; But actual x values are
     ; [..., -1.42, -0.42, 0.58, 1.58, ...]
     xarrc[*,i] -= xcc[i]
     yarrc[*,i] -= ycc[i]
     ; New grid coordinate system, recentered near the reference peak. In these
     ; coordinates, peak is at [XPEAK, YPEAK] - floor(XPEAK - 1, YPEAK - 1). So, e.g.,
     ; [0.42, 0.08] if [XPEAK, YPEAK] = [23.42, 28.08].
     ; But actual x values are:
     ; [-22, ..., -1, 0, 1, ...]
     xnewarr[*,i] -= floor(xcc[i])
     ynewarr[*,i] -= floor(ycc[i])

  endfor

  ; New x, y ranges. Pad with some additional columns / rows.
  newxran = double([round(min(xarrc))-1,round(max(xarrc))+1])
  newcols = fix(newxran[1]-newxran[0]+2)
  newyran = double([round(min(yarrc))-1,round(max(yarrc))+1])
  newrows = fix(newyran[1]-newyran[0]+2)

;; Enforce odd # of columns and rows
;  if not newcols then newcols+=1
;  if not newrows then newrows+=1
; Center galaxy in new grid ...
  if ~ keyword_set(nocenter) then begin
     coloff = round(double(newcols/2d))-1
     rowoff = round(double(newrows/2d))-1
; ... or not.
  endif else begin
;     coloff = fix(mean(xcc))+2
;     rowoff = fix(mean(ycc))+2
     coloff = -min(xnewarr)+2
     rowoff = -min(ynewarr)+2
  endelse
  ; this adds in some integer offset to the new coordinate system in every
  ; exposure.
  xnewarr += coloff
  ynewarr += rowoff
  crpix1 += coloff
  crpix2 += rowoff
  if keyword_set(crpix) then begin
     crpix1 = crpix[0]
     crpix2 = crpix[1]
  endif
  if keyword_set(crval) then begin
     crval1 = crval[0]
     crval2 = crval[1]
  endif

; Passing header structure to SXADDPAR won't update parameters. Something
; about how structure tags get passed to procedure and returned/updated
  pouthead=outheader.phu
  douthead=outheader.dat
  vouthead=outheader.var
  qouthead=outheader.dq

  if ~keyword_set(nophu) then begin
     sxaddpar,pouthead,'CRPIX1',crpix1
     sxaddpar,pouthead,'CRPIX2',crpix2
     if keyword_set(crval) then begin
        sxaddpar,pouthead,'CRVAL1',crval1
        sxaddpar,pouthead,'CRVAL2',crval2
     endif
  endif
  sxaddpar,douthead,'CRPIX1',crpix1
  sxaddpar,douthead,'CRPIX2',crpix2
  sxaddpar,vouthead,'CRPIX1',crpix1
  sxaddpar,vouthead,'CRPIX2',crpix2
  sxaddpar,qouthead,'CRPIX1',crpix1
  sxaddpar,qouthead,'CRPIX2',crpix2
  if keyword_set(crval) then begin
     sxaddpar,douthead,'CRVAL1',crval1
     sxaddpar,douthead,'CRVAL2',crval2
     sxaddpar,vouthead,'CRVAL1',crval1
     sxaddpar,vouthead,'CRVAL2',crval2
     sxaddpar,qouthead,'CRVAL1',crval1
     sxaddpar,qouthead,'CRVAL2',crval2
  endif
  
; Interpolating data to full grid and combining exposures.

; DATA

  print,'IFSR_MOSAIC: Interpolating data to full grid and combining.'

; Initialize arrays for new cube.
  datnewall = dblarr(newcols,newrows,nz0,nexp)-noval
  dqnewall = bytarr(newcols,newrows,nz0,nexp)+1b
  datnew = dblarr(newcols,newrows,nz0)
  dqnew = intarr(newcols,newrows,nz0)

  for i=0,nz0-1 do begin
     for j=0,nexp-1 do begin
        ; What we do here is index each original array such that the center of 
        ; each spaxel has a coordinate offset so that the peak reference gets shifted
        ; to the center of a spaxel upon interpolation.
        ; For xarr, when it indexes the original data array:
        ; [0, ..., 21, 22, 23, ...]
        ; But actual assigned x values are:
        ; [0.42, ..., 21.42, 22.42, 23.42, ...]
        ; in the new array, these correspond to
        ; [-22, ..., -1, 0, 1, ...] + some offset. So coordinate 22.42 gets
        ; mapped to coordinate 0 + some offset.
        ; The new coordinate in the final pixel array will have a unity offset, 
        ; of course.
        datnewall[xnewarr[*,j],ynewarr[*,j],i,j] = $
           interpolate(datall[*,*,i,j],$
                       xarr[*,j],yarr[*,j],$
                       /grid,missing=noval)
        tmp = interpolate(dqall[*,*,i,j],$
                          xarr[*,j],yarr[*,j],$
                          /grid,missing=0)
        zind = where(tmp le 0.01,ctz)
        oind = where(tmp gt 0.01,cto)
        if ctz gt 0 then tmp[zind] = 0
        if cto gt 0 then tmp[oind] = 1
        dqnewall[xnewarr[*,j],ynewarr[*,j],i,j] = tmp
     endfor
  endfor
  
  for i=0,newcols-1 do begin
     for j=0,newrows-1 do begin
        for k=0,nz0-1 do begin
           good = where(datnewall[i,j,k,*] ne noval AND $
                        dqnewall[i,j,k,*] eq 0,ctg)
           if ctg gt 0 then begin
              datnew[i,j,k] = median(datnewall[i,j,k,good],/double,/even)
              dqnew[i,j,k] = 0
           endif else begin
              datnew[i,j,k] = 0
              dqnew[i,j,k] = 1
           endelse
        endfor
     endfor
  endfor

  if ~keyword_set(nophu) then $
     mwrfits,cube.phu,outfile,pouthead,/create,silent=quiet
  mwrfits,datnew,outfile,douthead,create=~appenddat,silent=quiet

  datnewall=0
  datnew=0  

; VARIANCE

  print,'IFSR_MOSAIC: Interpolating variance to full grid and combining.'

  varnewall = dblarr(newcols,newrows,nz0,nexp)+noval
  varnew = dblarr(newcols,newrows,nz0)

  for i=0,nz0-1 do begin
     for j=0,nexp-1 do begin
        varnewall[xnewarr[*,j],ynewarr[*,j],i,j] = $
           interpolate(varall[*,*,i,j],$
                       xarr[*,j],yarr[*,j],$
                       /grid,missing=noval)
     endfor
  endfor
  for i=0,newcols-1 do begin
     for j=0,newrows-1 do begin
        for k=0,nz0-1 do begin
           good = where(varnewall[i,j,k,*] ne noval AND $
                        dqnewall[i,j,k,*] eq 0,ctg)
           if ctg gt 0 then begin
;              varnew[i,j,k] = median(varnewall[i,j,k,good],/double,/even)
              varnew[i,j,k] = total(varnewall[i,j,k,good],/double)/double(ctg)^2d
           endif else begin
              varnew[i,j,k] = 0
           endelse
        endfor
     endfor
  endfor
  

  mwrfits,varnew,outfile,vouthead,silent=quiet
  mwrfits,dqnew,outfile,qouthead,silent=quiet

  varnewall=0
  dqnewall=0
  varnew=0
  dqnew=0

end
