; docformat = 'rst'
;
;+
;
; Open a FITS data cube, call IFSR_MAKELINEMAP to create a linemap, and write the
; linemap to disk.
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
;      2018nov04, DSNR, created
;      2019jan25, DSNR, separated IFSR_MAKELINEMAP into subroutine
;      2019mar04, DSNR, moved IFSR_MAKELINEMAP to separate file
;
; :Copyright:
;    Copyright (C) 2018-19 David S. N. Rupke
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
pro ifsr_linemap,infile,outfile,wavesum,wavesub=wavesub,$
                 datext=datext,varext=varext,dqext=dqext

   bad=1d99

   if ~ keyword_set(datext) then datext=1
   if ~ keyword_set(varext) then varext=2
   if ~ keyword_set(dqext) then dqext=3

   header=1
   cube = ifsf_readcube(infile,header=header,datext=datext,varext=varext,$
                        dqext=dqext,/quiet)

   linesum = ifsr_makelinemap(cube,wavesum,wavesub=wavesub)

   appenddat=0b
   if datext eq 1 then begin
      writefits,outfile,[],header.phu
      appenddat=1b
   endif
   writefits,outfile,linesum.dat,header.dat,append=appenddat
   writefits,outfile,linesum.var,header.var,/append

end
