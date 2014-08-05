; docformat = 'rst'
;
;+ 
; 
; Merge data cubes processed with NIFS_ALIGN into a single cube using averaging,
; DQ rejection, and image quality thresholds from peak finding.
;
; :Categories:
;    IFSRED/NIFS
;
; :Returns:
;    Merged NIFS data cube.
;
; :Params:
;    inlist: in, required, type=string
;      List of input FITS files. Should be the same as outpsf for NIFS_ALIGN.
;    outfile: in, required, type=string
;      Output FITS file.
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
;      2013mar01, DSNR, copied from GMOS routine
;      2014aug05, DSNR, added documentation, license, and copyright
;
; :Copyright:
;    Copyright (C) 2013-2014 David S. N. Rupke
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
pro nifs_merge,inlist,outfile

;--------------------------------------------------------------------
; User inputs
;--------------------------------------------------------------------

  ncols = 45
  nrows = 45
  nz = 2040
; reference spectral coordinates
  crvalref = 18970.7d
  cdeltref = 2.135355d
  crpixref = 1
  wavelo = 19100
  wavehi = 23250
; maximum peak FWHM (") for inclusion in merge
  fwhmthresh = 0.14d

; compute reference spectral coordinates
  pix = dindgen(nz)
  waveref = crvalref + cdeltref*(pix+1-crpixref) 

; output array
  cube_dat_out = dblarr(ncols,nrows,nz)
  cube_var_out = dblarr(ncols,nrows,nz)
  cube_dq_out = dblarr(ncols,nrows,nz)

; Get filenames
  nfile=0
  openr,lun,inlist,/get_lun
  line=''
  while (~ EOF(lun)) do begin
     readf,lun,line
     linesplit = strsplit(line,/extract)
     avgfwhm = linesplit[5]
     if avgfwhm le fwhmthresh then begin
        if nfile eq 0 then filename=linesplit[1] $
        else filename=[filename,linesplit[1]]
        nfile++
     endif
  endwhile
  free_lun,lun
; Get headers
  dumy = mrdfits(filename[0],0,phu,/silent)
  dumyimg = mrdfits(filename[0],1,hdr,/silent)
; Cycle through wavelengths
  imodfiftylast=0
  for i=0,nz-1 do begin
     if waveref[i] ge wavelo AND waveref[i] le wavehi then begin
        slice_dat_tmp = dblarr(ncols,nrows,nfile)
        slice_var_tmp = dblarr(ncols,nrows,nfile)
        used_tmp = intarr(ncols,nrows,nfile)+1
        for j=0,nfile-1 do begin
           slice_dat_tmp_single = $
              mrdfits(filename[j],1,range=[i,i],/silent)
           slice_var_tmp_single = $
              mrdfits(filename[j],2,range=[i,i],/silent)
           slice_dq_tmp_single = $
              mrdfits(filename[j],3,range=[i,i],/silent)
           slice_dat_tmp[*,*,j] = slice_dat_tmp_single
           slice_var_tmp[*,*,j] = slice_var_tmp_single
;       Find zeroed out and bad regions
           used_tmp_single = intarr(ncols,nrows)+1 
           iunused_tmp_single = where(slice_dat_tmp_single eq 0)
           ;; iunused_tmp_single = where(slice_dat_tmp_single eq 0 OR $
           ;;                            slice_dq_tmp_single eq 1)
           used_tmp_single[iunused_tmp_single] = 0
           used_tmp[*,*,j] = used_tmp_single
        endfor
        slice_dat = total(slice_dat_tmp,3)/total(used_tmp,3)
        slice_var = total(slice_var_tmp,3)/total(used_tmp,3)/total(used_tmp,3)
        cube_dat_out[*,*,i] = slice_dat
        cube_var_out[*,*,i] = slice_var
     endif
     imodfifty = floor(double(i)/50d)
     if imodfifty ne imodfiftylast then print,'Slice ',i,' completed.'
     imodfiftylast = imodfifty
  endfor
  
  writefits,outfile,dumy,phu
  writefits,outfile,cube_dat_out,hdr,/append
  writefits,outfile,cube_var_out,hdr_var,/append
  ;; writefits,outfile,cube_dq_out,hdr_dq,/append
  
  cube_dat_out = 0

end
