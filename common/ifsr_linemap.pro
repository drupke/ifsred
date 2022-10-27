; docformat = 'rst'
;
;+
;
; Open a FITS data cube, call IFSR_MAKELINEMAP to create a linemap,
; and write the linemap to disk.
; 
; For now, subsampling does not correctly adjust the WCS of the input image.
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Output FITS file with linemap and variance.
;
; :Params:
;    infile: in, required, type=string
;      Path and filename of input file.
;    outfile: in, required, type=string
;      Path and filename of output file.
;    wavesum: in, required, type=dblarr(2)
;      Lower and upper limits of wavelength region over which to sum.
;
; :Keywords:
;    allowneg: in, optional, type=boolean
;      Default is to calculate continuum fluxes using positive values
;      only; this allows negative fluxes in the computation as well.
;    subsamplefac: in, optional, type=double
;      Factor by which to interpolate to a smaller grid.
;    subsampleout: out, optional, type=string
;      Paht and filename of output file if subsampling.
;    wavesub: in, optional, type=dblarr(2,N)
;      Lower and upper wavelength regions for subtraction continuum.
;      E.g., [5000,5100] or [[4800,4900],[5000,5100]].
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
;      2018nov04, DSNR, created
;      2019jan25, DSNR, separated IFSR_MAKELINEMAP into subroutine
;      2019mar04, DSNR, moved IFSR_MAKELINEMAP to separate file
;      2019may22, DSNR, keyword ALLOWNEG allows negative fluxes in computing 
;                       continuum region averages to subtract
;
; :Copyright:
;    Copyright (C) 2018--2019 David S. N. Rupke
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
pro ifsr_linemap,infile,outfile,wavesum,allowneg=allowneg,$
                 subsampleout=subsampleout,$
                 subsamplefac=subsamplefac,$
                 wavesub=wavesub,$
                 datext=datext,varext=varext,dqext=dqext


   bad=1d99

   if ~ keyword_set(datext) then datext=1
   if ~ keyword_set(varext) then varext=2
   if ~ keyword_set(dqext) then dqext=3

   header=1
   cube = ifsf_readcube(infile,header=header,datext=datext,varext=varext,$
                        dqext=dqext,/quiet)

   linesum = ifsr_makelinemap(cube,wavesum,wavesub=wavesub,allowneg=allowneg)

;  Subsampling interpolation
   if keyword_set(subsampleout) then begin
      if ~ keyword_set(subsamplefac) then subsamplefac=10d
      xresamp = (dindgen(cube.ncols*subsamplefac)-subsamplefac/2d + 0.5d)/subsamplefac
      yresamp = (dindgen(cube.nrows*subsamplefac)-subsamplefac/2d + 0.5d)/subsamplefac
      interp = interpolate(linesum.dat,xresamp,yresamp,/double,/grid,cubic=-0.5)
      writefits,subsampleout,interp,header.dat
   endif

   appenddat=0b
   if datext eq 1 then begin
      writefits,outfile,[],header.phu
      appenddat=1b
   endif
   writefits,outfile,linesum.dat,header.dat,append=appenddat
   writefits,outfile,linesum.var,header.var,/append
   

end
