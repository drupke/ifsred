; docformat = 'rst'
;
;+
;
; :Categories:
;    IFSRED/KCWI
;
; :Returns:
;
; :Params:
;    inraw: in, required, type=string
;    indat: in, required, type=string
;    inrep: in, required, type=string
;    minrepval: in, required, type=double
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
;      2018aug08  DSNR  created
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
pro ifsr_kcwisatcor,inraw,indat,inrep,minrepval,usebackup=usebackup,$
                    ampmode=ampmode,intscl=intscl
                    
;  ampmode
   if ~ keyword_set(ampmode) then ampmode = 9

;  Raw data
   rawdat = mrdfits(inraw+'.fits',0,rawhead,/silent,/unsigned)
;  Data sections
   rawsec1 = sxpar(rawhead,'DSEC1',/silent)
   rawsec2 = sxpar(rawhead,'DSEC2',/silent)
;  Convert data section strings to 4-element vectors, in IDL 0-offset array coords
   rawsec1p = kcwi_parse_imsec(rawsec1)-1
   rawsec2p = kcwi_parse_imsec(rawsec2)-1

;  new array; assume sections have same y coordinate range
   dx = (rawsec1p[1]-rawsec1p[0]+1) + (rawsec2p[0]-rawsec2p[1]+1)
   dy = rawsec1p[3]-rawsec1p[2]+1
   rawdatuse = uintarr(dx,dy)
;  endpoint of first half of new array, in new array coords
   x1end = rawsec1p[1]-rawsec1p[0]

;  populate new array
   rawdatuse[0:x1end,*] = rawdat[rawsec1p[0]:rawsec1p[1],rawsec1p[2]:rawsec1p[3]]
   rawdatuse[x1end+1:dx-1,*] = rawdat[rawsec2p[1]:rawsec2p[0],rawsec2p[2]:rawsec2p[3]]
   if ampmode eq 9 then rawdatuse = reverse(rawdatuse,2)

   ibad = where(rawdatuse ge minrepval,ctbad)
   iarrbad = array_indices(rawdatuse,ibad) ; For CR search (below)

   if ctbad gt 0 then begin

      if keyword_set(usebackup) then backstr='_nosatcor' else backstr = ''
      intdat = mrdfits(indat+'_int'+backstr+'.fits',0,inthead,/silent)
      mskdat = mrdfits(indat+'_msk'+backstr+'.fits',0,mskhead,/silent)
      vardat = mrdfits(indat+'_var'+backstr+'.fits',0,varhead,/silent)

      intrep = mrdfits(inrep+'_int.fits',0,rephead,/silent)
      mskrep = mrdfits(inrep+'_msk.fits',0,/silent)
      varrep = mrdfits(inrep+'_var.fits',0,/silent)

      intexp = sxpar(inthead,'XPOSURE',/silent)
      repexp = sxpar(rephead,'XPOSURE',/silent)
      if ~ keyword_set(intscl) then intscl = double(intexp)/double(repexp)

;     Look up and down CRPAD pixels (in columns) from each saturated point
;     to see if pixels near saturated point mistaken for CR. 
;     Replace these points as well.
      crpad = 4
      s = size(iarrbad)
      for i=0,ctbad-1 do begin
         if ((iarrbad[0,i] ge crpad)&&(iarrbad[0,i] le s[1]-crpad))>0 then begin
            mskdat_tmp = mskdat[iarrbad[0,i]-crpad:iarrbad[0,i]+crpad,iarrbad[1,i]]
            icr = where(mskdat_tmp eq 4,ctcr)

            if ctcr gt 0 then begin
               icr += ((iarrbad[1,i])*dx+iarrbad[0,i]-crpad)
               ibad = cgsetunion(ibad,icr)
            endif
         endif
      endfor

      intdat[ibad] = intrep[ibad]*intscl
      vardat[ibad] = varrep[ibad]*intscl
      mskdat[ibad] = mskrep[ibad]
      
      if ~ keyword_set(usebackup) then begin
         file_copy,indat+'_int.fits',indat+'_int_nosatcor.fits'
         file_copy,indat+'_msk.fits',indat+'_msk_nosatcor.fits'
         file_copy,indat+'_var.fits',indat+'_var_nosatcor.fits'
      endif

      writefits,indat+'_int.fits',intdat,inthead
      writefits,indat+'_msk.fits',mskdat,mskhead
      writefits,indat+'_var.fits',vardat,varhead
      
   endif else begin

      print,'No saturated pixels found.'

   endelse

end
