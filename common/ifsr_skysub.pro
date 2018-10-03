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
;    infile: in, required, type=string
;      Path and filename of input file.
;    outfile: in, required, type=string
;      Path and filename of output file.
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
;      2018sep18, DSNR, created
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
pro ifsr_skysub,infile,outfile,sumpar,nophu=nophu,spaxlist=spaxlist,$
                ignorepar=ignorepar

  header=1
  if keyword_set(nophu) then begin
     datext=-1
     varext=1
     dqext=2
     appenddat=0b
  endif else begin
     datext=1
     varext=2
     dqext=3
     appenddat=1b
  endelse
  cube = ifsf_readcube(infile,header=header,/quiet,datext=datext,varext=varext,$
                       dqext=dqext)  

;  Code copied from IFSR_SPAXSUM for defining / ignoring regions
   dx = cube.ncols
   dy = cube.nrows
   dz = cube.nz

   map_x = rebin(dindgen(dx)+1,dx,dy)
   map_y = rebin(transpose(dindgen(dy)+1),dx,dy)
   if keyword_set(spaxlist) then begin
      isum = []
      for i=0,n_elements(spaxlist[*,0])-1 do begin
         itmp = where(map_x eq spaxlist[i,0] AND $
                      map_y eq spaxlist[i,1])
         isum = [isum,itmp]
      endfor
   endif else begin
      if n_elements(sumpar) eq 3 then begin
;        Circular aperture
         map_r = sqrt((map_x - double(sumpar[0]))^2d + $
                      (map_y - double(sumpar[1]))^2d)
         isum = where(map_r le double(sumpar[2]))
      endif else if n_elements(sumpar) eq 4 then begin
;        Square aperture
         isumx = where(map_x ge double(sumpar[0]) AND $
                       map_x le double(sumpar[2]))
         isumy = where(map_y ge double(sumpar[1]) AND $
                       map_y le double(sumpar[3]))
         isum = cgsetintersection(isumx,isumy)
      end else message,'Number of elements in SUMPAR must be 3 or 4.'
   endelse

;  circular or square region to ignore
   if keyword_set(ignorepar) then begin
      if n_elements(ignorepar) eq 3 then begin
;        Circular aperture
         map_r = sqrt((map_x - double(ignorepar[0]))^2d + $
                      (map_y - double(ignorepar[1]))^2d)
         iignore = where(map_r le double(ignorepar[2]))
      endif else if n_elements(ignorepar) eq 4 then begin
;        Square aperture
         iignorex = where(map_x ge double(ignorepar[0]) AND $
                          map_x le double(ignorepar[2]))
         iignorey = where(map_y ge double(ignorepar[1]) AND $
                          map_y le double(ignorepar[3]))
         iignore = cgsetintersection(iignorex,iignorey)
      end else message,'Number of elements in IGNOREPAR must be 3 or 4.'
      isum = cgsetdifference(isum,iignore)
   endif
   
   skydatmed=dblarr(cube.nz)
   skydatavg=dblarr(cube.nz)
;   skyvar=dblarr(cube.nz)
   for i=0,cube.nz-1 do begin
      skydattmp = cube.dat[*,*,i]
;      skyvartmp = cube.var[*,*,i]
      skydqtmp = 1+(-1*cube.dq[*,*,i])
      skydattmp*=double(skydqtmp)
;      skyvartmp*=double(skydqtmp)
      nsky = n_elements(isum)
      ibad = where(skydqtmp eq 0,ctbad)
      if ctbad gt 0 then nsky-=ctbad
      skydatavg[i] = total(skydattmp[isum])/nsky
      skydatmed[i] = median(skydattmp[isum])
;      skyvar[i] = total(skyvartmp[isum])/nsky
   endfor      

   skydatarr = rebin(reform(skydatmed,1,1,dz),dx,dy,dz)
;   skyvararr = rebin(reform(skyvar,1,1,dz),dx,dy,dz)
  
   cube.dat -= skydatarr
;   cube.var += skyvararr

; Write output

   if ~ keyword_set(nophu) then writefits,outfile,cube.phu,header.phu
   writefits,outfile,cube.dat,header.dat,append=appenddat
   writefits,outfile,cube.var,header.var,/append
   writefits,outfile,cube.dq,header.dq,/append

end
