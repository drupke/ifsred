; docformat = 'rst'
;
;+
;
; Flip and/or rotate a data cube. Note that flips are performed before
; rotations.
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Flipped / rotated cube.
;
; :Params:
;    infile: in, required, type=string
;      Path and filename of input file.
;    outfile: in, required, type=string
;      Path and filename of output file.
;
;
; :Keywords:
;    xflip: in, optional, type=byte
;      Flip horizontally (left-to-right).
;    yflip: in, optional, type=byte
;      Flip vertically (top-to-bottom).
;    rotang: in, optional, type=double
;      Angle by which to rotate counterclockwise (must be 0, 90, 180,
;      or 270).
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
;      2014may23, DSNR, created
;
; :Copyright:
;    Copyright (C) 2014 David S. N. Rupke
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
pro ifsr_fliprot,infile,outfile,xflip=xflip,yflip=yflip,rotang=rotang

;  Read cube
   header=1
   cube = ifsf_readcube(infile,header=header,/quiet)  

   if keyword_set(xflip) then begin
      cube.dat = reverse(cube.dat,1)
      cube.var = reverse(cube.var,1)
      cube.dq = reverse(cube.dq,1)
   endif

   if keyword_set(yflip) then begin
      cube.dat = reverse(cube.dat,2)
      cube.var = reverse(cube.var,2)
      cube.dq = reverse(cube.dq,2)
   endif

   if keyword_set(rotang) then begin
      if rotang eq 0 then dir=0 $
      else if rotang eq 90 then dir=1 $
      else if rotang eq 180 then dir=2 $
      else if rotang eq 270 then dir=3 $
      else begin
         print,'IFSF_FLIPROT: Invalid parameter ROTANG.'
         exit
      endelse
      for i=0,cube.nz-1 do begin
         cube.dat[*,*,i] = rotate(cube.dat[*,*,i],dir)
         cube.var[*,*,i] = rotate(cube.var[*,*,i],dir)
         cube.dq[*,*,i] = rotate(cube.dq[*,*,i],dir)
      endfor
   endif

;  Write cube
   writefits,outfile,cube.phu,header.phu
   writefits,outfile,cube.dat,header.dat,/append
   writefits,outfile,cube.var,header.var,/append
   writefits,outfile,cube.dq,header.dq,/append

end
