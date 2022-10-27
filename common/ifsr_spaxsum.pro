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
;    outfile: in, required, type=string
;      Path and filename of output file.
;    sumpar: in, required, type=dblarr(3 or 4)
;      Parameters for summation, in unity-offset coordinates. Circular 
;      aperture: [x0, y0, r]. Rectangular aperture: [x1,y1,x2,y2]. Ignored if 
;      keyword SPAXLIST set.
;      
; :Keywords:
;    bunit: in, optional, type=str
;      Set data header keyword to BUNIT. Set var header keyword to (BUNIT)**2.
;    datext: in, optional, type=integer, default=1
;      Extension # of data plane. Set to a negative number if the correct
;      extension is 0, since an extension of 0 ignores the keyword.
;    dqext: in, optional, type=integer, default=3
;      Extension # of DQ plane. Set to a negative number if there is no DQ; DQ
;      plane is then set to 0.
;    exact: in, optional, type=boolean
;      If set and aperture is circular, compute and apply fractional pixel area
;      using exact circular region. Mirrors logic from ASTROLIB routine APER.PRO.
;      If set and a 2d spectral image, computes exact aperture including fractional
;      spaxels. Center of a spaxel is an integer.
;    fluxnorm: in, optional, type=double
;      Multiply data by this value.
;    fluxscale: in, optional, type=double
;      Divide data by this value. FLUXSCAL keyword added to header.
;    invvar: in, optional, type=boolean
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
;    varext: in, optional, type=integer, default=2
;      Extension # of variance plane.
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
;      2020may05, DSNR, added ability to sum 2D images; 
;         added CUNIT1 and BUNIT to outputs
;      2021dec07, DSNR, option to compute EXACT circular aperture
;      2022jan07, DSNR, now can compute exact circular aperture to ignore
;      2022jan21, DSNR, added _extra parameter to pass extra parameters 
;         to IFSF_READCUBE; added ability to compute exact aperture
;         for 2d spectral images
;      2022jan27, DSNR, added FLUXSCALE, FLUXNORM, BUNIT parameters; 
;        changed BUNIT in variance  to replicate BUNIT in variance header, 
;        not data header; fixed small bug in setting outheaddat for case of PHU
;      2022jul28, DSNR, removed NOPHU and added DATEXT, VAREXT, DQEXT; now
;        handles case of absent VAR, DQ planes
;
; :Copyright:
;    Copyright (C) 2015--2022 David S. N. Rupke
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
                 ignorepar=ignorepar,reversevardq=reversevardq,$
                 waveext=waveext,invvar=invvar,spaxarea=spaxarea,rotpa=rotpa,$
                 rotcent=rotcent,exact=exact,fluxscale=fluxscale,$
                 fluxnorm=fluxnorm,bunit=bunit,$
                 datext=datext,varext=varext,dqext=dqext,_extra=_extra
                 
   if ~ keyword_set(datext) then datext=1
   if ~ keyword_set(varext) then varext=2
   if ~ keyword_set(dqext) then dqext=3
   if keyword_set(reversevardq) then begin
      varext=3
      dqext=2
   endif

   if ~ keyword_set(invvar) then invvar=0
   header=1b
   if ~ keyword_set(waveext) then begin
      cube = ifsf_readcube(infile,/quiet,datext=datext,varext=varext,dqext=dqext,$
                           invvar=invvar,header=header,_extra=_extra)
   endif else begin
      cube = ifsf_readcube(infile,/quiet,datext=datext,varext=varext,dqext=dqext,$
                           waveext=waveext,invvar=invvar,header=header,_extra=_extra)
   endelse
   ; check if var and dq were found
   novar=0b
   if n_elements(cube.var) eq 1 then $
      if cube.var[0] eq -1 then $
         novar=1b
   nodq=0b
   if n_elements(cube.dq) eq 1 then $
      if cube.dq[0] eq -1 then $
         nodq=1b
   dx = cube.ncols
   dy = cube.nrows

   ; isum contains pixels to sum, either partially or completely
   ; iignore contains pixels to ignore, either partially or completely

   ; Sum in 2 dimensions
   if dy gt 1 then begin

      ; map of x- and y-coordinates, single-offset
      map_x = rebin(dindgen(dx) + 1d,dx,dy)
      map_y = rebin(transpose(dindgen(dy) + 1d),dx,dy)

      ; start masks with 0s
      ; fracmask is fractional flux of spaxel to include in aperture
      fracmask = dblarr(dx,dy)
      ; fracmask is fractional flux of spaxel to *not* include in aperture
      ignorefracmask = dblarr(dx,dy)

      ; List of specific spaxels to sum
      if keyword_set(spaxlist) then begin
         ; One spaxel
         if size(spaxlist,/n_dimensions) eq 1 then begin
            isum = where(map_x eq spaxlist[0] AND $
                         map_y eq spaxlist[1])
         ; >1 spaxel
         endif else begin
            isum = []
            for i=0,n_elements(spaxlist[*,0])-1 do begin
               itmp = where(map_x eq spaxlist[i,0] AND $
                            map_y eq spaxlist[i,1])
               isum = [isum,itmp]
            endfor
         endelse
         fracmask[isum] = 1d

      ; Circular or square aperture
      endif else begin

         ; Circular aperture
         if n_elements(sumpar) eq 3 then begin
            map_x_wrt_cent = abs(map_x - double(sumpar[0]))
            map_y_wrt_cent = abs(map_y - double(sumpar[1]))
            map_r = sqrt(map_x_wrt_cent^2d + map_y_wrt_cent^2d)

            ; Exact computation
            if keyword_set(exact) then begin

               ; Pixels fully in aperture, with fracmask = 1
               smallrad = sumpar[2]/sqrt(2d) - 0.5d
               igood = where(map_x_wrt_cent lt smallrad AND $
                  map_y_wrt_cent lt smallrad, ctgood)
               if ctgood GT 0 then fracmask[igood] = 1.0

               ; Pixels fully out of aperture, with fracmask = -1 (only used in
               ; next step)
               bigrad = sumpar[2] + 0.5d
               ibad = where(map_x_wrt_cent gt bigrad OR $
                  map_y_wrt_cent GT bigrad)
               fracmask[ibad] = -1

               ; Other pixels, for which to compute the fraction 0 < frac < 1
               ifrac = where(fracmask EQ 0d, ctfrac)
               if ctfrac GT 0 then fracmask[ifrac] = $
                  pixwt(double(sumpar[0]),double(sumpar[1]),double(sumpar[2]),$
                     map_x[ifrac],map_y[ifrac]) > 0d
               ; Set non-aperture pixels back to fracmask = 0
               fracmask[ibad] = 0

               isum = where(fracmask gt 0d,nsum)

            ; Approximation to nearest pixel
            endif else begin               

               isum = where(map_r le double(sumpar[2]),nsum)
               fracmask[isum] = 1d

            endelse

         ; Square aperture
         endif else if n_elements(sumpar) eq 4 then begin

            isumx = where(map_x ge double(sumpar[0]) AND $
                          map_x le double(sumpar[2]))
            isumy = where(map_y ge double(sumpar[1]) AND $
                          map_y le double(sumpar[3]))
            isum = cgsetintersection(isumx,isumy)
            fracmask[isum] = 1d 

         end else message,'Number of elements in SUMPAR must be 3 or 4.'

      endelse

   ; Sum in only 1 dimension
   endif else begin

      print,'IFSR_SPAXSUM: Assuming dispersion direction along rows.'
      map_x = dindgen(dx)
      fracmask = dblarr(dx)
      ; this allows to specify exact spaxel fractions to sum
      ; center of a spaxel is an integer, first spaxel spans 0.5 to 1.5
      ; 
      if keyword_set(exact) then begin
         ; convert to coordinate where first spaxel is 0.0 to 1.0
         sumpar -= 0.5d
         ; total number of spaxels in mask
         nsum = floor(sumpar[1]) - floor(sumpar[0])+ 1
         ; indices of spaxels to mask
         isum = indgen(nsum)+floor(sumpar[0])
         if nsum gt 1 then begin
            ; fraction of first spaxel in mask to mask
            fracmasklo = double(ceil(sumpar[0])) - double(sumpar[0])
            ; fraction of last spaxel in mask to mask
            fracmaskhi = double(sumpar[1]) - double(floor(sumpar[1]))
            fracmask[isum] = 1d
            fracmask[isum[0]] = fracmasklo
            fracmask[isum[nsum-1]] = fracmaskhi
         endif else begin
            fracmask[isum] = double(sumpar[1] - sumpar[0])
         endelse
      endif else begin
         isum = indgen(sumpar[1]-sumpar[0]+1)+sumpar[0]-1
         fracmask[isum] = 1d
      endelse

   endelse

   ;  circular or square region to ignore
   if keyword_set(ignorepar) then begin
      ; Circular aperture
      if n_elements(ignorepar) eq 3 then begin
         map_x_wrt_cent = abs(map_x - double(ignorepar[0]))
         map_y_wrt_cent = abs(map_y - double(ignorepar[1]))
         map_r = sqrt(map_x_wrt_cent^2d + map_y_wrt_cent^2d)
         if keyword_set(exact) then begin            
            
            ; Pixels fully in ignore aperture, with ignorefracmask = 1
            smallrad = ignorepar[2]/sqrt(2d) - 0.5d
            igood = where(map_x_wrt_cent lt smallrad AND $
               map_y_wrt_cent lt smallrad, ctgood)
            if ctgood GT 0 then ignorefracmask[igood] = 1.0

            ; Pixels fully out of ignore aperture, with ignorefracmask = 0
            bigrad = ignorepar[2] + 0.5d
            ibad = where(map_x_wrt_cent gt bigrad OR map_y_wrt_cent GT bigrad)
            ignorefracmask[ibad] = -1

            ; Other pixels, for which to compute the ignore fraction 0 < frac < 1
            ifrac = where(ignorefracmask EQ 0d, ctfrac)
            if ctfrac GT 0 then ignorefracmask[ifrac] = $
               pixwt(double(ignorepar[0]),double(ignorepar[1]),double(ignorepar[2]),$
                  map_x[ifrac], map_y[ifrac]) > 0d
            ignorefracmask[ibad] = 0

            ; subtract ignorefracmask from fracmask. 
            ; fractional pixels will be treated appropriately
            fracmask -= ignorefracmask
            ; case where ignorefracmask overcorrects fractional spaxel
            ineg = where(fracmask lt 0d,ctneg)
            if ctneg gt 0 then fracmask[ineg] = 0d

            ; recreate isum to remove completely ignored points
            isum = where(fracmask gt 0d)

         endif else begin               

            iignore = where(map_r le double(ignorepar[2]),nsum)
            isum = cgsetdifference(isum,iignore)
            fracmask[iignore] = 0d

         endelse

      endif else if n_elements(ignorepar) eq 4 then begin

;        Square aperture
         iignorex = where(map_x ge double(ignorepar[0]) AND $
                          map_x le double(ignorepar[2]))
         iignorey = where(map_y ge double(ignorepar[1]) AND $
                          map_y le double(ignorepar[3]))
         iignore = cgsetintersection(iignorex,iignorey)
         isum = cgsetdifference(isum,iignore)
         fracmask[iignore] = 0d

      end else message,'Number of elements in IGNOREPAR must be 3 or 4.'

   endif

;  plot image, fracmask, and aperture
;  white -> fracmask = 1
;  black -> fracmask = 0
;  grey -> 0 < fracmask < 1
   if cube.ndim eq 3 then begin
      set_plot,'x'
      cgdisplay,aspect=double(dy)/double(dx)
      cgimage,bytscl(fracmask,min=0,max=1),pos=[0.05,0.05,0.95,0.95],backg='Black'
      cgplot,[0],[0],/nodata,xran=[0.5,double(dx)+0.5],$
         pos=[0.05,0.05,0.95,0.95],yran=[0.5,double(dy)+0.5],/noerase,$
         axiscolor='White',ticklen=1d
      ;  aperture
      if n_elements(sumpar) eq 3 then $
         tvcircle,sumpar[2],sumpar[0],sumpar[1],color='Red',thick=2
   endif

   ; optional weights
   if ~ keyword_set(weights) then weights=dblarr(n_elements(isum))+1d
   
   outdat=dblarr(cube.nz)
   outvar=dblarr(cube.nz)
   outdq=bytarr(cube.nz)
   fluxmult=1d
   if keyword_set(spaxarea) then fluxmult = spaxarea
   for i=0,cube.nz-1 do begin
      if dy gt 1 then begin
         outdattmp = cube.dat[*,*,i]
         if ~novar then outvartmp = cube.var[*,*,i]
         if ~nodq then outdqtmp = cube.dq[*,*,i]
      endif else begin
;        This assumes dispersion is along rows
         outdattmp = cube.dat[i,*]
         if ~novar then outvartmp = cube.var[i,*]
         if ~nodq then outdqtmp = cube.dq[i,*]
      endelse
;     Documentation of ROT doesn't specify that rotation center
;     is in zero-offset coordinates but I tested it. Rotation can't be 
;     around a non-integer pixel center. The pivot keyword is unnecessary
;     if ROTCENT is the center of the image.
;     Do nearest-neighbor interpolation if it's a bad pixel mask.
      if keyword_set(rotpa) AND keyword_set(rotcent) then begin
         outdattmp = $
            rot(outdattmp,rotpa,1d,rotcent[0]-1,rotcent[1]-1,cubic=-0.5,/pivot)
         if ~novar then outvartmp = $
            rot(outvartmp,rotpa,1d,rotcent[0]-1,rotcent[1]-1,cubic=-0.5,/pivot)
         if ~nodq then outdqtmp = $
            rot(outdqtmp,rotpa,1d,rotcent[0]-1,rotcent[1]-1,/pivot)
      endif
      outdat[i] = total(outdattmp[isum]*weights*fracmask[isum])*fluxmult
; I *think* dat/err should not change when going from S.B. to flux . ..
      if ~novar then $
         outvar[i] = total(outvartmp[isum]*weights*fracmask[isum])*fluxmult^2d
      if ~nodq then begin
         outdq[i] = total(outdqtmp[isum]*weights)
         if outdq[i] gt 0b then outdq[i] = 1b
      endif
   endfor      

;   fxhmake,outhead,/extend,/date
   
   if keyword_set(fluxnorm) then begin
      outdat *= fluxnorm
      if ~novar then outvar *= fluxnorm*fluxnorm
   endif
   if keyword_set(fluxscale) then begin
      outdat /= fluxscale
      if ~novar then outvar /= fluxscale*fluxscale
   endif
   
   if datext eq -1 then $
      fxhmake,outheaddat,outdat,/extend,/date $
   else $
      fxhmake,outheaddat,outdat,/xtension,/date
   sxaddpar,outheaddat,'EXTNAME','SCI'
   sxaddpar,outheaddat,'DISPAXIS',1
   sxaddpar,outheaddat,'WCSDIM',1
   sxaddpar,outheaddat,'CTYPE1','LINEAR'
   sxaddpar,outheaddat,'CRPIX1',cube.crpix
   sxaddpar,outheaddat,'CRVAL1',cube.crval
   sxaddpar,outheaddat,'CD1_1',cube.cdelt
   sxaddpar,outheaddat,'CDELT1',cube.cdelt
   if cube.cunit ne '' then sxaddpar,outheaddat,'CUNIT1',cube.cunit
   if keyword_set(bunit) then sxaddpar,outheaddat,'BUNIT',bunit $
   else if cube.bunit ne '' then sxaddpar,outheaddat,'BUNIT',cube.bunit
   if keyword_set(fluxscale) then sxaddpar,outheaddat,'FLUXSCAL',fluxscale
   
   if ~novar then begin
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
      if keyword_set(bunit) then sxaddpar,outheadvar,'BUNIT','('+bunit+')**2' $
      else if cube.bunit_var ne '' then sxaddpar,outheadvar,'BUNIT_VAR',cube.bunit
      if keyword_set(fluxscale) then sxaddpar,outheadvar,'FLUXSCAL',fluxscale*fluxscale
   endif

   if ~nodq then begin
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
   endif

   ; create PHU if needed
   if datext ne -1 then begin
      mwrfits,[],outfile,header.phu,/create
      create=0b
   endif else create=1b
   mwrfits,outdat,outfile,outheaddat,create=create
   if keyword_set(reversevardq) then begin
      if ~nodq then mwrfits,outdq,outfile,outheaddq
      if ~novar then mwrfits,outvar,outfile,outheadvar
   endif else begin
      if ~novar then mwrfits,outvar,outfile,outheadvar
      if ~nodq then mwrfits,outdq,outfile,outheaddq
   endelse
   if keyword_set(waveext) then $
      mwrfits,cube.wave,outfile,header_wave

end
