; docformat = 'rst'
;
;+
;
; :Categories:
;    IFSRED
;
; :Returns:
;
; :Params:
;    incube: in, required, type=string
;      Path and filename for input data cube. Cube should contain data, 
;      variance, and DQ planes in extensions 1-3.
;    outim: in, required, type
;    w1: in, required, type=integer
;      Degree of polynomial fit to centroid columns vs. wavelength.
;    w2: in, required, type=integer
;      Degree of polynomial fit to centroid rows vs. wavelength.
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
;      2017jan17, DSNR, removed invocation of SYMBOLS routine
;                       
; :Copyright:
;    Copyright (C) 2017 David S. N. Rupke
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
pro ifsr_stack,infile,outfile,w1,w2,mean=mean

;  Load data cube
   header=1
   cube = ifsf_readcube(infile,header=header,/quiet)

   iw1 = value_locate(cube.wave,w1)
   iw2 = value_locate(cube.wave,w2)
   outim = median(cube.dat[*,*,iw1:iw2],/double,dim=3)
   
   writefits,outfile,outim
   
end
