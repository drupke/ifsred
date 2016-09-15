; docformat = 'rst'
;
;+
;
;  Change the native OSIRIS data cube storage (wavelength, row, col) to match
;  the typical IFS data cube storage (col, row, wavlength). Note that this
;  preserves increasing column and row numbers as they were in the original
;  data cube. It also
;   -- adjusts the keywords accordingly
;   -- converts error into variance
;   -- sets good = 0 and bad = 1 in the DQ plane
;   -- makes data and variance data types double instead of float, and dq
;      data type int instead of unsigned int (the latter prevents IRAF issues
;      later and matches the GMOS DQ data type)
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Rearrange 
;
; :Params:
;    infile: in, required, type=str
;    outfile: in, required, type=str
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
;      2016sep12, DSNR, created
;
; :Copyright:
;    Copyright (C) 2016 David S. N. Rupke
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
pro osiris_transpose,infile,outfile

   dat = mrdfits(infile,0,header_dat,/silent)
   err = mrdfits(infile,1,header_err,/silent)
   dq = mrdfits(infile,2,header_err,/silent)

   dat_new = double(transpose(dat,[2,1,0]))
   var = double(transpose(err,[2,1,0]))^2d
   dq_new = fix(transpose(dq,[2,1,0]))

   igd = where(dq_new gt 0,ctgd)
   ibd = where(dq_new eq 0,ctbd)
   if ctgd gt 0 then dq_new[igd] = 0b
   if ctbd gt 0 then dq_new[ibd] = 1b

   header_dat_new=header_dat
   sxaddpar,header_dat_new,'CTYPE1',sxpar(header_dat,'CTYPE3',/silent)
   sxaddpar,header_dat_new,'CTYPE3',sxpar(header_dat,'CTYPE1',/silent)
   sxaddpar,header_dat_new,'CUNIT1',sxpar(header_dat,'CUNIT3',/silent)
   sxaddpar,header_dat_new,'CUNIT3',sxpar(header_dat,'CUNIT1',/silent)
   sxaddpar,header_dat_new,'CRVAL1',sxpar(header_dat,'CRVAL3',/silent)
   sxaddpar,header_dat_new,'CRVAL3',sxpar(header_dat,'CRVAL1',/silent)
   sxaddpar,header_dat_new,'CRPIX1',sxpar(header_dat,'CRPIX3',/silent)
   sxaddpar,header_dat_new,'CRPIX3',sxpar(header_dat,'CRPIX1',/silent)
   sxaddpar,header_dat_new,'CDELT1',sxpar(header_dat,'CDELT3',/silent)
   sxaddpar,header_dat_new,'CDELT3',sxpar(header_dat,'CDELT1',/silent)
;  Remove these orientation keywords, rather than adjust; not sure
;  how they're affected by the transposing.
   sxdelpar,header_dat_new,['PC1_1','PC2_2','PC2_3','PC3_2','PC3_3']

   writefits,outfile,dat_new,header_dat_new,/check
   writefits,outfile,var,header_err,/append,/check
   writefits,outfile,dq_new,header_dq,/append,/check
   
end
