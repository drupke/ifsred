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
; :Categories:
;    IFSRED
;
; :Returns:
;    Combined data cube. Header of cube is the same as the header of
;    the first exposure in the file list.
;    
; :Params:
;    infiles: in, required, type=strarr
;      Filenames for data cubes of individual exposures.
;    outfile: in, required, type=string
;      Path and filename of output file.
;
; :Keywords:
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
;
; :Copyright:
;    Copyright (C) 2014--2016 David S. N. Rupke
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
pro ifsr_mosaic,infiles,outfile,indir=indir,nocenter=nocenter,nophu=nophu

; Value for data outside the bounds of interpolation
  noval=-999d

  if not keyword_set(indir) then indir=''

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
     cube = ifsf_readcube(indir+infiles[i],header=header,/quiet,$
                          datext=datext,varext=varext,dqext=dqext)

;    Set output header to header of first cube
     if i eq 0 then outheader=header

     if ~ keyword_set(nophu) then begin
        writefits,outfile,cube.phu,outheader.phu
        pkheader = header.phu
        appenddat=1
     endif else begin
        pkheader = header.dat
        appenddat=0
     endelse
     
;    Fitted galaxy center
     xcc[i] = sxpar(pkheader,'XPEAK',/silent,count=ctx) - 1
     ycc[i] = sxpar(pkheader,'YPEAK',/silent,count=cty) - 1
     if ~ keyword_set(nophu) AND (ctx eq 0 OR cty eq 0) then begin
;       Check science header if it's missing from PHU
        xcc[i] = sxpar(header.dat,'XPEAK','exposure'+string(i+1,format='(I0)'),$
                       /silent) - 1
        ycc[i] = sxpar(header.dat,'YPEAK','exposure'+string(i+1,format='(I0)'),$
                       /silent) - 1
     endif
;    Residual from nearest integer pixel
     xccd[i] = xcc[i] - floor(xcc[i])
     yccd[i] = ycc[i] - floor(ycc[i])

;    Initialize arrays to hold all exposures, as well as some other 
;    intermediate arrays
     if i eq 0 then begin
        nz0 = cube.nz
        crpix0 = cube.crpix
        datall = dblarr(cube.ncols,cube.nrows,cube.nz,nexp)
        varall = dblarr(cube.ncols,cube.nrows,cube.nz,nexp)
        dqall = dblarr(cube.ncols,cube.nrows,cube.nz,nexp)
        xarr = rebin(dindgen(cube.ncols),cube.ncols,nexp)
        yarr = rebin(dindgen(cube.nrows),cube.nrows,nexp)
        xarrc = rebin(dindgen(cube.ncols),cube.ncols,nexp)
        yarrc = rebin(dindgen(cube.nrows),cube.nrows,nexp)
        xnewarr = rebin(indgen(cube.ncols),cube.ncols,nexp)
        ynewarr = rebin(indgen(cube.nrows),cube.nrows,nexp)
     endif

;    Put everything on a common wavelength solution, if it is not already.
     dcrpix = cube.crpix - crpix0
     if cube.nz le nz0 then begin
        if dcrpix le 0 then begin
           datall[*,*,abs(dcrpix):cube.nz-1,i] = $
              cube.dat[*,*,0:cube.nz-1-abs(dcrpix)]
           varall[*,*,abs(dcrpix):cube.nz-1,i] = $
              cube.var[*,*,0:cube.nz-1-abs(dcrpix)]
           dqall[*,*,abs(dcrpix):cube.nz-1,i] = $
              cube.dq[*,*,0:cube.nz-1-abs(dcrpix)]
        endif else begin
           datall[*,*,0:cube.nz-1-dcrpix,i] = cube.dat[*,*,dcrpix:cube.nz-1]
           varall[*,*,0:cube.nz-1-dcrpix,i] = cube.var[*,*,dcrpix:cube.nz-1]
           dqall[*,*,0:cube.nz-1-dcrpix,i] = cube.dq[*,*,dcrpix:cube.nz-1]
        endelse
     endif else begin
        if dcrpix le 0 then begin
           datall[*,*,abs(dcrpix):nz0-1,i] = cube.dat[*,*,0:nz0-1-abs(dcrpix)]
           varall[*,*,abs(dcrpix):nz0-1,i] = cube.var[*,*,0:nz0-1-abs(dcrpix)]
           dqall[*,*,abs(dcrpix):nz0-1,i] = cube.dq[*,*,0:nz0-1-abs(dcrpix)]
        endif else begin
           datall[*,*,0:nz0-1-dcrpix,i] = cube.dat[*,*,dcrpix:nz0-1]
           varall[*,*,0:nz0-1-dcrpix,i] = cube.var[*,*,dcrpix:nz0-1]
           dqall[*,*,0:nz0-1-dcrpix,i] = cube.dq[*,*,dcrpix:nz0-1]
        endelse
     endelse

;    Add in integer residual w.r.t. center
     xarr[*,i] += xccd[i]
     yarr[*,i] += yccd[i]
     xarrc[*,i] -= xcc[i]
     yarrc[*,i] -= ycc[i]
     xnewarr[*,i] -= floor(xcc[i])
     ynewarr[*,i] -= floor(ycc[i])

  endfor

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
     coloff = fix(mean(xcc))+2
     rowoff = fix(mean(ycc))+2
  endelse
  xnewarr += coloff
  ynewarr += rowoff
  

; Interpolating data to full grid and combining exposures.

; DATA

  print,'IFSR_MOSAIC: Interpolating data to full grid and combining.'

; Initialize arrays for new cube.
  datnewall = dblarr(newcols,newrows,cube.nz,nexp)-noval
  dqnewall = bytarr(newcols,newrows,cube.nz,nexp)+1b
  datnew = dblarr(newcols,newrows,cube.nz)
  dqnew = intarr(newcols,newrows,cube.nz)

  for i=0,cube.nz-1 do begin
     for j=0,nexp-1 do begin
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
        for k=0,cube.nz-1 do begin
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

  writefits,outfile,datnew,outheader.dat,append=appenddat

  datnewall=0
  datnew=0  

; VARIANCE

  print,'IFSR_MOSAIC: Interpolating variance to full grid and combining.'

  varnewall = dblarr(newcols,newrows,cube.nz,nexp)+noval
  varnew = dblarr(newcols,newrows,cube.nz)

  for i=0,cube.nz-1 do begin
     for j=0,nexp-1 do begin
        varnewall[xnewarr[*,j],ynewarr[*,j],i,j] = $
           interpolate(varall[*,*,i,j],$
                       xarr[*,j],yarr[*,j],$
                       /grid,missing=noval)
     endfor
  endfor
  for i=0,newcols-1 do begin
     for j=0,newrows-1 do begin
        for k=0,cube.nz-1 do begin
           good = where(varnewall[i,j,k,*] ne noval AND $
                        dqnewall[i,j,k,*] eq 0,ctg)
           if ctg gt 0 then begin
              varnew[i,j,k] = median(varnewall[i,j,k,good],/double,/even)
           endif else begin
              varnew[i,j,k] = 0
           endelse
        endfor
     endfor
  endfor
  
  writefits,outfile,varnew,outheader.var,/append
  writefits,outfile,dqnew,outheader.dq,/append

  varnewall=0
  dqnewall=0
  varnew=0
  dqnew=0

end
