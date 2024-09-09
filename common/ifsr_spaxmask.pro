; docformat = 'rst'
;
;+
;
; Mask spaxels for summing in a data cube over a square or circular aperture.
;
; :Categories:
;    IFSRED
 ;
; :Returns:
;    Rebinned data cube.
;
; :Params:
;    sumpar: in, required, type=dblarr(3 or 4)
;      Parameters for summation, in unity-offset coordinates. Circular 
;      aperture: [x0, y0, r]. Rectangular aperture: [x1,y1,x2,y2]. Ignored if 
;      keyword SPAXLIST set.
;      
; :Keywords:
;    exact: in, optional, type=boolean
;      If set and aperture is circular, compute and apply fractional pixel area
;      using exact circular region. Mirrors logic from ASTROLIB routine APER.PRO.
;      If set and a 2d spectral image, computes exact aperture including fractional
;      spaxels. Center of a spaxel is an integer.
;    spaxlist: in, optional, type=dblarr(Nspax,2)
;      List of spaxels to combine. If this is chosen, SUMPAR parameter is ignored.
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
;      2023mar31, DSNR, created
;
; :Copyright:
;    Copyright (C) 2023 David S. N. Rupke
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
function ifsr_spaxmask, dx, dy, sumpar, spaxlist=spaxlist,$
   ignorepar=ignorepar, exact=exact
                 

   ; isum contains pixels to sum in 2d, either partially or completely
   ; iignore contains pixels to ignore, either partially or completely

   ; Sum in 2 dimensions
   if dy gt 1 then begin

      ; map of x- and y-coordinates, single-offset
      map_x = rebin(dindgen(dx) + 1d,dx,dy)
      map_y = rebin(transpose(dindgen(dy) + 1d),dx,dy)

      ; start masks with 0s
      ; fracmask is fractional flux of spaxel to include in aperture
      fracmask = dblarr(dx,dy)
      ; ignorefracmask is fractional flux of spaxel to *not* include in aperture
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

      print,'IFSR_SPAXMASK: Assuming dispersion direction along rows.'
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

   return, {fracmask: fracmask, isum: isum}
   end