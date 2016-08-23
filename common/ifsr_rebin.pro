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
;    infile: in, required, type=string
;      Path and filename of input file.
;    outfile: in, required, type=string
;      Path and filename of output file.
;    factor: in, required, type=integer
;      Factor by which data cube is rebinned; in either dimension, 
;      old dimension / new dimension (not counting pixels dropped to make
;      rebinning occur by an integer number)
;    xlim: in, required, type=intarr(2)
;      Indices for lower and upper column numbers to rebin. First col and row
;      is (1,1).
;    ylim: in, required, type=intarr(2)
;      Indices for lower and upper row numbers to rebin. First col and row
;      is (1,1).
;      
; :Keywords:
;    applydq: in, option, type=byte
;      Use rebinned DQ plane to zero data and variance in output spaxels where all
;      input spaxels were bad.
;    fluxperarea: in, option, type=byte
;      Conserve flux per unit area, rather than total flux.
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
;      2014feb07, DSNR, created
;      2016feb03, DSNR, added way to apply DQ to output data and variance
;      2016jul13, DSNR, added option for conserving flux/area; previously
;                       routine conserved total flux, not flux/area
;
; :Copyright:
;    Copyright (C) 2014--2016 David S. N. Rupke
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
pro ifsr_rebin,infile,outfile,factor,xlim,ylim,applydq=applydq,$
               fluxperarea=fluxperarea

   factor = double(factor)

   header=1
   cube = ifsf_readcube(infile,header=header)

   nx = fix((xlim[1]-xlim[0]+1) / factor)
   ny = fix((ylim[1]-ylim[0]+1) / factor)

;  Switch to 0-offset coords
   xlim--
   ylim--

;  Rebin. This command does not re-scale flux; so, e.g., rebinning a 2x2
;  array with 1s to a 1x1 array returns 1.
; 
;  The fluxfactor should be factor^2 if total flux is being conserved.
;  E.g., rebinning a 2x2 array with 1s to a 1x1 array should return 0.25.
;  If total flux per unit area is being conserved, then rebinning a 2x2 array 
;  with 1s to a 1x1 array should return 1. Thus, when flux per unit area is
;  multiplied by area, then the correct total flux is returned in either case.
   if keyword_set(fluxperarea) then fluxfact=1d else fluxfact = factor^2d
   dat_rb = rebin(cube.dat[xlim[0]:xlim[1],ylim[0]:ylim[1],*],$
                  nx,ny,cube.nz) * fluxfact
   var_rb = rebin(cube.var[xlim[0]:xlim[1],ylim[0]:ylim[1],*],$
                  nx,ny,cube.nz) * fluxfact
;  DQ is treated differently b/c not propagating flux, but a flag
   dq_rb  = fix(rebin(double(cube.dq[xlim[0]:xlim[1],ylim[0]:ylim[1],*]),$
                      nx,ny,cube.nz) * factor^2d)
;  Procedure is to flag an output spaxel as bad if one of the spaxels that went
;  into it is bad. However, also zero out data and variance if every one of spaxels
;  that went into binned spaxel were bad.
   if keyword_set(applydq) then begin
      fac_arr = dat_rb*0d + factor^2d/(factor^2d - double(dq_rb))
      ibad = where(dq_rb eq factor^2d)
      fac_arr[ibad] = 0d
      dat_rb *= fac_arr
      var_rb *= fac_arr
   endif
   bad = where(dq_rb gt 0,ct)
   dq_rb[bad] = 1
;  dq_rb -= 32768

   newheader_phu = header.phu
   newheader_dat = header.dat
   newheader_var = header.var
   newheader_dq = header.dq
   cd11 = sxpar(newheader_phu,'CD1_1',/silent)
   cd12 = sxpar(newheader_phu,'CD1_2',/silent)
   cd21 = sxpar(newheader_phu,'CD2_1',/silent)
   cd22 = sxpar(newheader_phu,'CD2_2',/silent)
   sxaddpar,newheader_phu,'CD1_1',factor*cd11
   sxaddpar,newheader_phu,'CD1_2',factor*cd12
   sxaddpar,newheader_phu,'CD2_1',factor*cd21
   sxaddpar,newheader_phu,'CD2_2',factor*cd22
   sxaddpar,newheader_dat,'CD1_1',-factor/10d
   sxaddpar,newheader_dat,'CD2_2',-factor/10d
   sxaddpar,newheader_dat,'CDELT1',-factor/10d
   sxaddpar,newheader_dat,'CDELT2',-factor/10d
   sxaddpar,newheader_var,'CD1_1',-factor/10d
   sxaddpar,newheader_var,'CD2_2',-factor/10d
   sxaddpar,newheader_var,'CDELT1',-factor/10d
   sxaddpar,newheader_var,'CDELT2',-factor/10d
   sxaddpar,newheader_dq,'CD1_1',-factor/10d
   sxaddpar,newheader_dq,'CD2_2',-factor/10d
   sxaddpar,newheader_dq,'CDELT1',-factor/10d
   sxaddpar,newheader_dq,'CDELT2',-factor/10d
 
;  Write output

   writefits,outfile,cube.phu,newheader_phu
   writefits,outfile,dat_rb,newheader_dat,/append
   writefits,outfile,var_rb,newheader_var,/append
   writefits,outfile,dq_rb,newheader_dq,/append

end
