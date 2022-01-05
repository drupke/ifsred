; docformat = 'rst'
;
;+
;
; Aperture photometry for an IFS datacube.
; 
; :Categories:
;    IFSRED
;
; :Returns:
;   
;
; :Params:
;    filename: in, required, type=string
;    apr: in, required, type=dblarr(N)
;    lamrange: in, required, type=dblarr(2)
;      Lower and upper limits of wavelength regions to include in fit.
;    
; :Keywords:
;    cent: in, optional, type=dblarr(2)
;      Center in X and Y; single-offset coordinates.
;    quiet: in, optional, type=byte
;      Suppress output of centroids and half-widths to terminal.
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
;      2021nov24, DSNR, added NaN rejection
;                       
; :Copyright:
;    Copyright (C) 2021 David S. N. Rupke
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
function ifsr_aper,filename,apr,lamrange,cent=cent,quiet=quiet,datext=datext,$
   varext=varext,dqext=dqext,funit=funit

   if ~keyword_set(datext) then datext=1
   if ~keyword_set(varext) then varext=2
   if ~keyword_set(dqext) then dqext=3
   if ~keyword_set(funit) then funit=1d

   header=1b
   cube = ifsf_readcube(filename,quiet=quiet,datext=datext,varext=varext,$
      dqext=dqext,header=header)

   if not keyword_set(cent) then begin
      cent=dblarr(2)
      if datext eq -1 then $
         cent = [sxpar(header.dat,'XPEAK'),sxpar(header.dat,'YPEAK')] $
      else $
         cent = [sxpar(header.phu,'XPEAK'),sxpar(header.phu,'YPEAK')]
   endif

;  Indices of fit range in wavelength
   indrange = intarr(2)
   indrange[0] = value_locate(cube.wave,lamrange[0])
   indrange[1] = value_locate(cube.wave,lamrange[1])

   wavsub = cube.wave[indrange[0]:indrange[1]]
   datsub = cube.dat[*,*,indrange[0]:indrange[1]]
   varsub = cube.var[*,*,indrange[0]:indrange[1]]
   dqsub = cube.dq[*,*,indrange[0]:indrange[1]]
   nlam = n_elements(wavsub)

;  Zero out DQed or NANed data
   ibad = where(dqsub eq 1 OR ~ finite(datsub))
   datsub[ibad] = 0d
   varsub[ibad] = 0d

;  Sum over wavelength space
   dlam = wavsub[1:nlam-1]-wavsub[0:nlam-2]
   dlam = [dlam[0],dlam]
   dlamcube = rebin(reform(dlam,1,1,nlam),cube.ncols,cube.nrows,nlam)
   datsum = total(datsub*dlamcube,3,/double)*funit
   varsum = total(varsub*dlamcube^2,3,/double)*funit^2

   ;  Normalize to the same # of wavelength pixels summed at each spaxel
   fluxct = (indrange[1]-indrange[0]+1) - total(dqsub,3,/double)
   datsum /= fluxct/max(fluxct)
   varsum /= fluxct/max(fluxct)

;  Flag spaxels with no data
   inoflux = where(fluxct eq 0,ctnoflux)
   if ctnoflux gt 0 then begin
     datsum[inoflux]=0d
     varsum[inoflux]=1d99
   endif

   aper,datsum,[cent[0]],[cent[1]],fluxes,dumy0,dumy1,dumy2,1d,apr,[0d,0d],[0d,0d],$
      setskyval=0d,/flux,/silent

   return,fluxes

end
