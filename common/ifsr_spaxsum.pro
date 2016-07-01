; docformat = 'rst'
;
;+
;
; Sum spaxels in a data cube over a square or circular aperture. Center must be 
; in integer coordinates, and for a square aperture the full widths must be
; even numbers.
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
;      Parameters for summation. First two elements are center for summation, 
;      in unity-offset x and y spaxel coordinates. If circular aperture,
;      third element is radius in spaxels and there is no fourth element. If
;      square aperture, then third and fourth elements are aperture full widths
;      in spaxels in x and y directions. Ignored if keyword SPAXLIST set.
;      
; :Keywords:
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
;      2015jan17, DSNR, copied from GMOS_NUCSPEC_MRK231
;      2015may15, DSNR, added SPAXLIST option
;
; :Copyright:
;    Copyright (C) 2015 David S. N. Rupke
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
pro ifsr_spaxsum,infile,outfile,sumpar,spaxlist=spaxlist,weights=weights

   cube = ifsf_readcube(infile,/quiet)
   dx = cube.ncols
   dy = cube.nrows

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
         map_r = sqrt((map_x - sumpar[0])^2d + (map_y - sumpar[1])^2d)
         isum = where(map_r le sumpar[2])
      endif else if n_elements(sumpar) eq 4 then begin
;        Square aperture
         isumx = where(map_x ge sumpar[0]-sumpar[2]/2 AND $
                       map_x le sumpar[0]+sumpar[2]/2)
         isumy = where(map_y ge sumpar[0]-sumpar[3]/2 AND $
                       map_y le sumpar[0]+sumpar[3]/2)
         isum = cgsetintersection(isumx,isumy)
      end else begin
         print,'IFSR_SPAXSUM: Number of elements in SUMPAR must be 3 or 4.'
         exit
      endelse
   endelse
   
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

   fxhmake,outhead,/extend,/date

   fxhmake,outheaddat,outdat,/xtension,/date
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

   writefits,outfile,[],outhead
   writefits,outfile,outdat,outheaddat,/append
   writefits,outfile,outvar,outheadvar,/append
   writefits,outfile,outdq,outheaddq,/append

end
