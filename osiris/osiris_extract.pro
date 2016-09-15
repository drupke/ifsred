; docformat = 'rst'
;
;+
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Fits file for creating a single spectrum from a square region in OSIRIS data.
;
; :Params:
;    infile: in, required, type=str
;    outfile: in, required, type=str
;    xrange: in, required, type=dblarr(2)
;      X range of spaxels to sum.
;    yrange: in, required, type=dblarr(2)
;      Y range of spaxels to sum.
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
;      2013apr29, DSNR, created
;      2016sep08, DSNR, documented, copyrighted, licensed
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
pro osiris_extract,infile,outfile,xrange,yrange

  dat = mrdfits(infile,0,dat_hdr,/silent)
  err = mrdfits(infile,1,err_hdr,/silent) 
  dq = mrdfits(infile,2,dq_hdr,/silent)

  var = err^2d
  dat_out = total(dat[*,*,xrange[0]:xrange[1]],3)
  dat_out = total(dat_out[*,yrange[0]:yrange[1]],2)
  var_out = total(var[*,*,xrange[0]:xrange[1]],3)
  var_out = total(var_out[*,yrange[0]:yrange[1]],2)
  dq_out = total(dq[*,*,xrange[0]:xrange[1]],3)
  dq_out = total(dq_out[*,yrange[0]:yrange[1]],2)
  err_out = sqrt(var_out)
  dq_out /= dq_out
  dq_out *= 9
  dq_out = byte(dq_out)

  writefits, outfile, dat_out, dat_hdr
  writefits, outfile, err_out, /append
  writefits, outfile, dq_out, /append

end
