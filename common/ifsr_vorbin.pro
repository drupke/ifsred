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
;      2018dec29, DSNR, created
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
pro ifsr_vorbin,invorfile,incube,outcube,$
                datext=datext,varext=varext

   bad=1d99
   if ~ keyword_set(datext) then datext=1
   if ~ keyword_set(varext) then varext=2

   restore,file=invorfile

   header=1 
   cube = ifsf_readcube(incube,header=header,/quiet,$
                        datext=datext,varext=varext,dqext=-1)  
   newcube = cube.dat*0d
   newvar = cube.var*0d

   binarr = indgen(cube.ncols,cube.nrows) + !VALUES.D_NaN
   xind = vordat.xin + vordat.center[0] - 1d
   yind = vordat.yin + vordat.center[1] - 1d
   nbins = max(vordat.binNum)+1
   for i=0,nbins-1 do begin
      ibin = where(vordat.binNum eq i,ctbin)
      if ctbin eq 1 then begin
         binarr[xind[ibin[0]],yind[ibin[0]]] = vordat.binNum[ibin[0]] + 1
         newcube[xind[ibin[0]],yind[ibin[0]],*] = $
            cube.dat[xind[ibin[0]],yind[ibin[0]],*]
         newvar[xind[ibin[0]],yind[ibin[0]],*] = $
            cube.var[xind[ibin[0]],yind[ibin[0]],*]
      endif else begin
         bindat = dblarr(cube.nz)
         binvar = dblarr(cube.nz)
         for j=0,ctbin-1 do begin
            binarr[xind[ibin[j]],yind[ibin[j]]] = vordat.binNum[ibin[j]] + 1
            bindat += cube.dat[xind[ibin[j]],yind[ibin[j]],*]
            binvar += cube.var[xind[ibin[j]],yind[ibin[j]],*]
         endfor
         bindat /= double(ctbin)
         binvar /= double(ctbin)
         for j=0,ctbin-1 do begin
            newcube[xind[ibin[j]],yind[ibin[j]],*] = bindat
            newvar[xind[ibin[j]],yind[ibin[j]],*] = binvar
         endfor
      endelse
   endfor

;  Write output
   create=1b
   if datext eq 1 then begin
      mwrfits,cube.phu,outcube,header.phu,/create
      create=0b
   endif
   mwrfits,newcube,outcube,header.dat,create=create
   mwrfits,newvar,outcube,header.var
   mwrfits,binarr,outcube

end
