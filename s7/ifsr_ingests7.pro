; docformat = 'rst'
;
;+
;
; :Categories:
;    IFSRED/S7
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
;      2018apr14  DSNR  created
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
pro ifsr_ingests7,inblue,inred,outfile,voronoi=voronoi,scalevar=scalevar

   bheader=1b
   rheader=1b
   cubeblue = ifsf_readcube(inblue,header=bheader,datext=-1,varext=1,dqext=-1)
   cubered = ifsf_readcube(inred,header=rheader,datext=-1,varext=1,dqext=-1)
   vormap = mrdfits(inblue,'VORONOI_BIN_MAP',maphead)

;  Get new wavelength array. Set up rebinning on red side.
   stitchwave = 5600d
   newdelt = double(cubeblue.cdelt)
;  index of stitch wavelength on blue side
   iwavehiblue = value_locate(cubeblue.wave,stitchwave)
;  new wavelength array on blue side
   newwave = cubeblue.wave[0:iwavehiblue]
;  first wavelength of new red array
   newwavelored = cubeblue.wave[iwavehiblue]+newdelt
;  last wavelength of old red array
   oldwavehired = cubered.wave[cubered.nz-1]
;  approx. width of new array
   tmpwaveredrange = double(oldwavehired - newwavelored)
;  number of elements in new array
   nnewwavered = floor(tmpwaveredrange / newdelt)+1d
;  exact width of new red array
   newwaveredrange = double(nnewwavered-1)*newdelt
;  new array
   newwavered = newwavelored + dindgen(nnewwavered)/$
                (double(nnewwavered)-1d)*newwaveredrange
   newwave = [newwave,newwavered]
    
;;  Get new wavelength array. Set up rebinning on both sides.
;   stitchwave = 5000d
;;  new dispersion is twice old red dispersion
;   newdelt = double(cubered.cdelt)*2d
;   wavelored = value_locate(cubered.wave,stitchwave)
;   iwavehired = cubered.nz-1
;;  ensure there are an even number of points on the red side for binning
;   noldwavered = fix(iwavehired-iwavelored+1)
;   if noldwavered then iwavehired--
;   nnewwavered = iwavehired - iwavelored + 1
;   newwavered = cubered.wave[iwavelored:iwavehired]
;   newwavehiblue = cubered.wave[iwavelored]-newdelt
;   oldwaveloblue = cubeblue.wave[0]
;   tmpwavebluerange = double(newwavehiblue - oldwaveloblue)
;   nnewwaveblue = floor(tmpwavebluerange / newdelt)+1
;   newwaveloblue = newwavehiblue - newdelt*(double(nnewwaveblue)-1d)
;   newwavebluerange = newwavehiblue - newwaveloblue
;   newwaveblue = newwaveloblue + dindgen(nnewwaveblue)/$
;                 (double(nnewwaveblue)-1d)*newwavebluerange
;   newwave = [newwaveblue,newwavered]
;   newnz = nnewwaveblue+nnewwavered
;
;;  bin up red data
;   newdatred = rebin(cubered.dat[*,*,iwavelored:iwavehired],
;   newvarred = dblarr(cubeblue.ncols,cubeblue.nrows,nnewwave)

;  resample red data
   if keyword_set(voronoi) then begin
      
      newdatred = dblarr(cubeblue.ncols,nnewwavered)
      newvarred = dblarr(cubeblue.ncols,nnewwavered)

      for i=0,cubered.ncols-1 do begin
         newdatred[i,*] = $
            interpol(cubered.dat[i,*],cubered.wave,newwavered,$
                     /nan,/spline)
         newvarred[i,*] = $
            interpol(cubered.var[i,*],cubered.wave,newwavered,$
                     /nan,/spline)
      endfor

;     adjust variance downward by sampling
      newvarred*=cubered.cdelt/newdelt
      
      newdat = [[cubeblue.dat[*,0:iwavehiblue]],[newdatred]]
      newvar = [[cubeblue.var[*,0:iwavehiblue]],[newvarred]]
    
      newheader_dat = bheader.dat
      newheader_var = bheader.var

      sxaddpar,newheader_dat,'CRVAL2',cubeblue.wave[0]
      sxaddpar,newheader_var,'CRVAL2',cubeblue.wave[0]
      sxaddpar,newheader_dat,'CRPIX2',1
      sxaddpar,newheader_var,'CRPIX2',1
      sxaddpar,newheader_dat,'CDELT2',newdelt
      sxaddpar,newheader_var,'CDELT2',newdelt
      sxaddpar,newheader_dat,'C2_2',newdelt
      sxaddpar,newheader_var,'C2_2',newdelt

   endif else begin
      
      newdatred = dblarr(cubeblue.ncols,cubeblue.nrows,nnewwavered)
      newvarred = dblarr(cubeblue.ncols,cubeblue.nrows,nnewwavered)

      for i=0,cubered.ncols-1 do begin
         for j=0,cubered.nrows-1 do begin
            newdatred[i,j,*] = $
               interpol(cubered.dat[i,j,*],cubered.wave,newwavered,$
                        /nan,/spline)
            newvarred[i,j,*] = $
               interpol(cubered.var[i,j,*],cubered.wave,newwavered,$
                        /nan,/spline)
         endfor
      endfor
        
      newdat = [[[cubeblue.dat[*,*,0:iwavehiblue]]],[[newdatred]]]
      newvar = [[[cubeblue.var[*,*,0:iwavehiblue]]],[[newvarred]]]
    
      newheader_dat = bheader.dat
      newheader_var = bheader.var

      sxaddpar,newheader_dat,'CRVAL3',cubeblue.wave[0]
      sxaddpar,newheader_var,'CRVAL3',cubeblue.wave[0]
      sxaddpar,newheader_dat,'CRPIX3',1
      sxaddpar,newheader_var,'CRPIX3',1
      sxaddpar,newheader_dat,'CDELT3',newdelt
      sxaddpar,newheader_var,'CDELT3',newdelt
      sxaddpar,newheader_dat,'C3_3',newdelt
      sxaddpar,newheader_var,'C3_3',newdelt
      
   endelse

   if keyword_set(scalevar) then newvar *= scalevar

   writefits,outfile,newdat/1d-15,newheader_dat
   writefits,outfile,newvar/1d-30,newheader_var,/append
   writefits,outfile,vormap,maphead,/append

end
