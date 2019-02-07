; docformat = 'rst'
;
;+
;
; Sum spaxels in a data cube over a square or circular aperture.
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
;    nophu: in, optional, type=byte
;      Added flag to indicate that the 0th extension contains the data, not
;      the PHU.
;    outfile: in, required, type=string
;      Path and filename of output file.
;    sumpar: in, required, type=dblarr(3 or 4)
;      Parameters for summation, in unity-offset coordinates. Circular 
;      aperture: [x0, y0, r]. Rectangular aperture: [x1,y1,x2,y2]. Ignored if 
;      keyword SPAXLIST set.
;      
; :Keywords:
;    fluxmult: in, optional, type=double
;      Multiplier for flux in each spaxel.
;    invvar: in, optional, type=byte
;      Set if the data cube holds inverse variance instead of variance. The
;      output structure will still contain the variance.
;    spaxlist: in, optional, type=dblarr(Nspax,2)
;      List of spaxels to combine. If this is chosen, SUMPAR parameter is ignored.
;    waveext: in, optional, type=integer
;      The extention number of a wavelength array.
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
;      2015jan17, DSNR, copied from GMOS_NUCSPEC_MRK231
;      2015may15, DSNR, added SPAXLIST option
;      2016sep12, DSNR, added NOPHU option
;      2018feb08, DSNR, added WAVEEXT and INVVAR keywords
;
; :Copyright:
;    Copyright (C) 2015--2018 David S. N. Rupke
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
pro ifsr_spaxsum,infile,outfile,sumpar,spaxlist=spaxlist,weights=weights,$
                 nophu=nophu,ignorepar=ignorepar,reversevardq=reversevardq,$
                 waveext=waveext,invvar=invvar
                 
   if keyword_set(nophu) then begin
      datext=-1
      varext=1
      dqext=2
      if keyword_set(reversevardq) then begin
         varext=2
         dqext=1
      endif
      appenddat=0b
   endif else begin
      datext=1
      varext=2
      dqext=3
      if keyword_set(reversevardq) then begin
         varext=3
         dqext=2
      endif
      appenddat=1b
   endelse
   if ~ keyword_set(waveext) then waveext=0
   if ~ keyword_set(invvar) then invvar=0
   cube = ifsf_readcube(infile,/quiet,datext=datext,varext=varext,dqext=dqext,$
                        waveext=waveext,invvar=invvar)
   dx = cube.ncols
   dy = cube.nrows

   map_x = rebin(dindgen(dx) + 1d,dx,dy)
   map_y = rebin(transpose(dindgen(dy) + 1d),dx,dy)
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
   
   if ~ keyword_set(weights) then weights=dblarr(n_elements(isum))+1d
   
   outdat=dblarr(cube.nz)
   outvar=dblarr(cube.nz)
   outdq=bytarr(cube.nz)
   for i=0,cube.nz-1 do begin
      outdattmp = cube.dat[*,*,i]
      outvartmp = cube.var[*,*,i]
      outdqtmp = cube.dq[*,*,i]
      outdat[i] = total(outdattmp[isum]*weights)
      outvar[i] = total(outvartmp[isum]*weights)
      outdq[i] = total(outdqtmp[isum]*weights)
      if outdq[i] gt 0b then outdq[i] = 1b
   endfor      

;   fxhmake,outhead,/extend,/date

   if ~ keyword_set(nophu) then $
      fxhmake,outheaddat,outdat,/xtension,/date $
   else $
      fxhmake,outheaddat,outdat,/extend,/date
   sxaddpar,outheaddat,'EXTNAME','SCI'
   sxaddpar,outheaddat,'DISPAXIS',1
   sxaddpar,outheaddat,'WCSDIM',1
   sxaddpar,outheaddat,'CTYPE1','LINEAR'
   sxaddpar,outheaddat,'CRPIX1',cube.crpix
   sxaddpar,outheaddat,'CRVAL1',cube.crval
   sxaddpar,outheaddat,'CD1_1',cube.cdelt
   sxaddpar,outheaddat,'CDELT1',cube.cdelt

   fxhmake,outheadvar,outvar,/xtension,/date
   sxaddpar,outheadvar,'EXTNAME','VAR'
   sxaddpar,outheadvar,'DISPAXIS',1
   sxaddpar,outheadvar,'WCSDIM',1
   sxaddpar,outheadvar,'CTYPE1','LINEAR'
   sxaddpar,outheadvar,'CRPIX1',cube.crpix
   sxaddpar,outheadvar,'CRVAL1',cube.crval
   sxaddpar,outheadvar,'CD1_1',cube.cdelt
   sxaddpar,outheadvar,'CDELT1',cube.cdelt

   fxhmake,outheaddq,outdq,/xtension,/date
   sxaddpar,outheaddq,'EXTNAME','SCI'
   sxaddpar,outheaddq,'DISPAXIS',1
   sxaddpar,outheaddq,'WCSDIM',1
   sxaddpar,outheaddq,'CTYPE1','LINEAR'
   sxaddpar,outheaddq,'CRPIX1',cube.crpix
   sxaddpar,outheaddq,'CRVAL1',cube.crval
   sxaddpar,outheaddq,'CD1_1',cube.cdelt
   sxaddpar,outheaddq,'CDELT1',cube.cdelt

   if ~ keyword_set(nophu) then writefits,outfile,[],outhead
   writefits,outfile,outdat,outheaddat,append=appenddat
   if keyword_set(reversevardq) then begin
      writefits,outfile,outdq,outheaddq,/append
      writefits,outfile,outvar,outheadvar,/append
   endif else begin
      writefits,outfile,outvar,outheadvar,/append
      writefits,outfile,outdq,outheaddq,/append
   endelse
   if keyword_set(waveext) then $
      writefits,outfile,cube.wave,header_wave,/append

end
