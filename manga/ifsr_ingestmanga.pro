; docformat = 'rst'
;
;+
;
;
; :Categories:
;    IFSRED/MANGA
;
; :Returns:
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
;      2018feb28  DSNR  created
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
;    along with this program.  If not, see http://www.gnu.org/licenses/.
;
;-
pro ifsr_ingestmanga,infile,outfile

    header=1b
    cube = ifsf_readcube(infile,header=header,waveext=4,gooddq=[0,32,64,96],$
                         /invvar,/linearize)

;   For some reason SXADDPAR doesn't work when applied directly to structures, 
;   so create new variables.
    newheader_dat = header.dat
    newheader_var = header.var
    newheader_dq = header.dq

    sxaddpar,newheader_dat,'CRVAL3',cube.crval
    sxaddpar,newheader_var,'CRVAL3',cube.crval
    sxaddpar,newheader_dq,'CRVAL3',cube.crval
    sxaddpar,newheader_dat,'CRPIX3',cube.crpix
    sxaddpar,newheader_var,'CRPIX3',cube.crpix
    sxaddpar,newheader_dq,'CRPIX3',cube.crpix
    sxaddpar,newheader_dat,'CDELT3',cube.cdelt
    sxaddpar,newheader_var,'CDELT3',cube.cdelt
    sxaddpar,newheader_dq,'CDELT3',cube.cdelt
    sxaddpar,newheader_dat,'C3_3',cube.cdelt
    sxaddpar,newheader_var,'C3_3',cube.cdelt
    sxaddpar,newheader_dq,'C3_3',cube.cdelt
    writefits,outfile,cube.phu,header.phu
    writefits,outfile,cube.dat,newheader_dat,/append
    writefits,outfile,cube.var,newheader_var,/append
    writefits,outfile,cube.dq,newheader_dq,/append

end
