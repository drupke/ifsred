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
;    invvar: in, optional, type=byte
;      Set if the data cube holds inverse variance instead of variance. The
;      output structure will still contain the variance.
;    rotpa: in, optional, type=double
;      PA to rotate image before extracting. For, e.g., a long-slit aperture,
;      select PA of aperture in degrees E of N. Image will rotate in other direction
;      but extraction will be in the correct sense.
;    rotcent: in, optional, type=dblarr(2)
;      Center of rotation, in integer single-offset coordinates
;    spaxarea: in, optional, type=double
;      If the input fluxes are per arcsec^2, this keyword specifies the
;      area in square arcseconds for a spaxel so that the output is in flux
;      rather than surface brightness.
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
;      2019jun28, DSNR, added SPAXAREA keyword
;      2019jul08, DSNR, added ROTPA and ROTCENT keywords
;      2020may05, DSNR, added ability to sum 2D images; added CUNIT1 and BUNIT to outputs
;
; :Copyright:
;    Copyright (C) 2015--2020 David S. N. Rupke
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
                 waveext=waveext,invvar=invvar,spaxarea=spaxarea,rotpa=rotpa,$
                 rotcent=rotcent
                 
   if keyword_set(nophu) then begin
      datext=-1
      varext=1
      dqext=2
      if keyword_set(reversevardq) then begin
         varext=2
         dqext=1
      endif
      create=1b
   endif else begin
      datext=1
      varext=2
      dqext=3
      if keyword_set(reversevardq) then begin
         varext=3
         dqext=2
      endif
      create=0b
   endelse

   if ~ keyword_set(invvar) then invvar=0
   header=1b
   if ~ keyword_set(waveext) then begin
      cube = ifsf_readcube(infile,/quiet,datext=datext,varext=varext,dqext=dqext,$
                           invvar=invvar,header=header)
   endif else begin
      cube = ifsf_readcube(infile,/quiet,datext=datext,varext=varext,dqext=dqext,$
                           waveext=waveext,invvar=invvar,header=header)
   endelse
   dx = cube.ncols
   dy = cube.nrows

   if dy gt 1 then begin
      map_x = rebin(dindgen(dx) + 1d,dx,dy)
      map_y = rebin(transpose(dindgen(dy) + 1d),dx,dy)
      if keyword_set(spaxlist) then begin
         if size(spaxlist,/n_dimensions) eq 1 then begin
            isum = where(map_x eq spaxlist[0] AND $
                         map_y eq spaxlist[1])
         endif else begin
            isum = []
            for i=0,n_elements(spaxlist[*,0])-1 do begin
               itmp = where(map_x eq spaxlist[i,0] AND $
                            map_y eq spaxlist[i,1])
               isum = [isum,itmp]
            endfor
         endelse
      endif else begin
         if n_elements(sumpar) eq 3 then begin
;           Circular aperture
            map_r = sqrt((map_x - double(sumpar[0]))^2d + $
                         (map_y - double(sumpar[1]))^2d)
            isum = where(map_r le double(sumpar[2]))
         endif else if n_elements(sumpar) eq 4 then begin
;           Square aperture
            isumx = where(map_x ge double(sumpar[0]) AND $
                          map_x le double(sumpar[2]))
            isumy = where(map_y ge double(sumpar[1]) AND $
                          map_y le double(sumpar[3]))
            isum = cgsetintersection(isumx,isumy)
         end else message,'Number of elements in SUMPAR must be 3 or 4.'
      endelse
   endif else begin
      print,'IFSR_SPAXSUM: Assuming dispersion direction along rows.'
      map_x = dindgen(dx)
      isum = dindgen(sumpar[1]-sumpar[0]+1)+sumpar[0]-1
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
   fluxmult=1d
   if keyword_set(spaxarea) then fluxmult = spaxarea
   for i=0,cube.nz-1 do begin
      if dy gt 1 then begin
         outdattmp = cube.dat[*,*,i]
         outvartmp = cube.var[*,*,i]
         outdqtmp = cube.dq[*,*,i]
      endif else begin
;        This assumes dispersion is along rows
         outdattmp = cube.dat[i,*]
         outvartmp = cube.var[i,*]
         outdqtmp = cube.dq[i,*]
      endelse
;     Documentation of ROT doesn't specify that rotation center
;     is in zero-offset coordinates but I tested it. Rotation can't be 
;     around a non-integer pixel center. The pivot keyword is unnecessary
;     if ROTCENT is the center of the image.
;     Do nearest-neighbor interpolation if it's a bad pixel mask.
      if keyword_set(rotpa) AND keyword_set(rotcent) then begin
         outdattmp = rot(outdattmp,rotpa,1d,rotcent[0]-1,rotcent[1]-1,cubic=-0.5,$
                         /pivot)
         outvartmp = rot(outvartmp,rotpa,1d,rotcent[0]-1,rotcent[1]-1,cubic=-0.5,$
                         /pivot)
         outdqtmp = rot(outdqtmp,rotpa,1d,rotcent[0]-1,rotcent[1]-1,$
                         /pivot)
      endif
      outdat[i] = total(outdattmp[isum]*weights)*fluxmult
; I *think* dat/err should not change when going from S.B. to flux . ..
      outvar[i] = total(outvartmp[isum]*weights)*fluxmult^2d
      outdq[i] = total(outdqtmp[isum]*weights)
      if outdq[i] gt 0b then outdq[i] = 1b
   endfor      

;   fxhmake,outhead,/extend,/date
   
   if ~ keyword_set(nophu) then $
      fxhmake,outhead,outdat,/xtension,/date $
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
   if cube.cunit ne '' then sxaddpar,outheaddat,'CUNIT1',cube.cunit
   if cube.bunit ne '' then sxaddpar,outheaddat,'BUNIT',cube.bunit

   fxhmake,outheadvar,outvar,/xtension,/date
   sxaddpar,outheadvar,'EXTNAME','VAR'
   sxaddpar,outheadvar,'DISPAXIS',1
   sxaddpar,outheadvar,'WCSDIM',1
   sxaddpar,outheadvar,'CTYPE1','LINEAR'
   sxaddpar,outheadvar,'CRPIX1',cube.crpix
   sxaddpar,outheadvar,'CRVAL1',cube.crval
   sxaddpar,outheadvar,'CD1_1',cube.cdelt
   sxaddpar,outheadvar,'CDELT1',cube.cdelt
   if cube.cunit ne '' then sxaddpar,outheadvar,'CUNIT1',cube.cunit
   if cube.bunit ne '' then sxaddpar,outheadvar,'BUNIT',cube.bunit

   fxhmake,outheaddq,outdq,/xtension,/date
   sxaddpar,outheaddq,'EXTNAME','SCI'
   sxaddpar,outheaddq,'DISPAXIS',1
   sxaddpar,outheaddq,'WCSDIM',1
   sxaddpar,outheaddq,'CTYPE1','LINEAR'
   sxaddpar,outheaddq,'CRPIX1',cube.crpix
   sxaddpar,outheaddq,'CRVAL1',cube.crval
   sxaddpar,outheaddq,'CD1_1',cube.cdelt
   sxaddpar,outheaddq,'CDELT1',cube.cdelt
   if cube.cunit ne '' then sxaddpar,outheaddq,'CUNIT1',cube.cunit
   if cube.bunit ne '' then sxaddpar,outheaddq,'BUNIT',cube.bunit

   if ~ keyword_set(nophu) then mwrfits,[],outfile,header.phu,/create
   mwrfits,outdat,outfile,outheaddat,create=create
   if keyword_set(reversevardq) then begin
      mwrfits,outdq,outfile,outheaddq
      mwrfits,outvar,outfile,outheadvar
   endif else begin
      mwrfits,outvar,outfile,outheadvar
      mwrfits,outdq,outfile,outheaddq
   endelse
   if keyword_set(waveext) then $
      mwrfits,cube.wave,outfile,header_wave

end
