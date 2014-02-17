function readspec,infile,extension=extension,header=header
;
; History
; 09april  DSNR  created (because there is nothing like it in Astrolib)
; 09jul01  DSNR  modified to include log-linear coordinates
; 11jan21  DSNR  added option to output header
; 13mar05  DSNR  modified treatment of error messages; not sure how
;                this will work ...
;

  errmsg=''
  
  if ~ keyword_set(extension) then extension=0

  flux = mrdfits(infile,extension,header,/silent)

  s = size(flux)
  if s[0] eq 1 then naps=1 else naps=s[2]
  npix = n_elements(flux)/naps
  wave = dindgen(npix)

  ct = intarr(6)
  ctype = sxpar(header,'CTYPE1',/silent,count=count)
  ct[0] = count
  crval = double(sxpar(header,'CRVAL1',/silent,count=count))
  ct[1] = count
  crpix = double(sxpar(header,'CRPIX1',/silent,count=count))
  ct[2] = count
  cdelt = double(sxpar(header,'CDELT1',/silent,count=count))
  ct[3] = count
  cd11  = double(sxpar(header,'CD1_1',/silent,count=count))
  ct[4] = count
  dcflag= sxpar(header,'DC-FLAG',/silent,count=count)
  ct[5] = count

  disp = cdelt
  nopar = where(ct eq 0,ctpar)
  if ctpar then begin
     if ctpar gt 1 then begin
     ;; if ctpar gt 1 OR (nopar[0] ne 3 AND nopar[0] ne 4) then begin
        errmsg = 'Header keywords missing.'
        goto,error
     endif else begin
        if nopar[0] eq 3 then disp = cd11
     endelse
  endif

  if (ctype eq 'LINEAR  ' OR ctype eq 'lambda  ') then begin
     wave = crval + disp*(wave-crpix+1)
     if dcflag then wave = 10d^wave
  endif else begin
     errmsg = 'Dispersion solution not linear.'
     goto,error
  endelse

error:
  if errmsg ne '' then print,'READSPEC: ERROR: ',errmsg,$
                             format='(T5,A,A)'

  return,[[wave],[flux]]

end
