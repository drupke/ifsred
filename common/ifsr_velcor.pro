; docformat = 'rst'
;
;+
;
; Convert from air to vacuum wavelengths and/or apply heliocentric correction.
; Spectrum linearized on return. Only good to nearest whole JD, 
;
; :Categories:
;    IFSRED
;
; :Returns:
;    Corrected spectrum or spectra.
;
; :Params:
;    infile: in, required, type=string
;      Path and filename of input file.
;    outfile: in, required, type=string
;      Path and filename of output file.
;      
; :Keywords:
;    vaccor: in, optional, type=boolean
;      If set, compute vacuum correction using AIRTOVAC from IDLAstro library.
;    heliocor: in, optional, type=boolean
;      If set, compute heliocentric correction using JULDAY and X_KECKHELIO. The
;      latter uses BARYVEL from IDLAstro.
;    cdelt: in, optional, type=double
;      Output linear dispersion. If not set, default to average value determined
;      from corrected wavelength array.
;    decimal: in, optional, type=boolean
;      Toggle if RA/DEC are input in decimal degrees.
;    mdyhms: in, optional, type=dblarr(3-6)
;      [Month, day, and year], and optionally [[[hour], minute], second] of 
;      observation input to JULDAY. Required if HELIOCOR set.
;    obs: in, optional, type=string, default='keck'
;      Observatory string from OBSERVATORY routine. 
;    ra: in, optional, type=string
;    decl: in, optional, type=string
;      Coordinates of target, in degrees or as a string sexagesimal. If not
;      set, looks for RA header keyword.
;    datext: in, optional, type=integer, default=1
;      Extension # of data plane. Set to a negative number if the correct
;      extension is 0, since an extension of 0 ignores the keyword.
;    varext: in, optional, type=integer, default=2
;      Extension # of variance plane.
;    dqext: in, optional, type=integer, default=3
;      Extension # of DQ plane. Set to a negative number if there is no DQ.
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
;      2022jan04, DSNR, created
;
; :Copyright:
;    Copyright (C) 2021 David S. N. Rupke
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
pro ifsr_velcor,infile,outfile,vaccor=vaccor,heliocor=heliocor,cdelt=cdelt,$
   decimal=decimal,mdyhms=mdyhms,obs=obs,ra=ra,decl=decl,$
   datext=datext,varext=varext,dqext=dqext

   if ~keyword_set(datext) then datext=1
   if ~keyword_set(varext) then varext=2
   if ~keyword_set(dqext) then dqext=3
   
   header=1
   cube = ifsf_readcube(infile,header=header,datext=datext,varext=varext,$
      dqext=dqext)

   if keyword_set(vaccor) OR keyword_set(heliocor) then begin
      wave = cube.wave
      dat = cube.dat
      var = cube.var
      dq = cube.dq 

      ; air to vacuum correction
      if keyword_set(vaccor) then begin
         print,'IFSR_VELCOR: Performing air-to-vaccum conversion.'
         airtovac,wave,newwave
      endif else newwave = wave

      ; heliocentric correction
      if keyword_set(heliocor) then begin
         if not keyword_set(obs) then begin
            print,'IFSR_VELCOR: Setting observatory location to Keck'
            obs = 'keck'
         endif
         ; Get RA and dec
         ; Assume header keywords are 'ra' and 'dec' in sexagesimal
         if not keyword_set(ra) AND not keyword_set(decl) then begin
            if datext eq 1 then useheader = header.phu $
            else useheader = header.dat
            radeg = sxpar(useheader,'ra')
            decdeg = sxpar(useheader,'dec')
         endif else begin
            radeg = ra
            decdeg = decl
         endelse
         ; if in decimal degrees
         if keyword_set(decimal) then begin
            radeg = double(radeg)
            decdeg = double(decdeg)
         endif else begin
            radeg = tenv(radeg)*15d
            decdeg = tenv(decdeg)
         endelse
         ; compute correction
         if not keyword_set(mdyhms) then begin
            print,'IFSR_VELCOR: Must specify date as [month, day, year].'
            print,'Optionally, add [[[hour], minute], second].'
            stop
         endif else begin
            hms = dblarr(3)
            for i=0,n_elements(mdyhms)-4 do hms[i] = mdyhms[3+i]
            jd = julday(mdyhms[0],mdyhms[1],mdyhms[2],hms[0],hms[1],hms[2])
            print,'IFSR_VELCOR: Julian day is ',string(jd,format='(D0.2)')
         endelse
         ; this routine produces the negative of what needs to be added to
         ; observed wavelength: dlam = lam * v_Earth,Sun/c,
         ; or lam_obs' = lam_obs + lam_obs * v_ES/c = lam_obs (1 + v_ES/c)
         rvcor = x_keckhelio(radeg,decdeg,jd=jd,obs=obs)
         print,'IFSR_VELCOR: Earth motion w.r.t. Sun/Earth barycenter is ',$
            string(-rvcor,format='(D0.2)'),' km/s in direction of target.'
         newwave *= 1d - rvcor/299792d
         print,'IFSR_VELCOR: Average wavelength change is ',$
            string(mean(newwave-wave),format='(D0.2)'),' A.'
      endif
      
      wave = newwave

      ; Re-linearize wavelengths
      ; average new cdelt
      if not keyword_set(cdelt) then $
         cdelt = double(wave[cube.nz-1]-wave[0])/double(cube.nz-1)
      ; non-linear wavelengths
      waveold = wave
      ; create new array with linear dispersion
      wave = double(wave[0]) + dindgen(cube.nz)*cdelt
      crpix = 1
      crval = wave[0]
      IF cube.ndim eq 3 then begin
         for i=0,cube.ncols-1 do begin
            for j=0,cube.nrows-1 do begin
               dat[i,j,*] = interpol(cube.dat[i,j,*],waveold,wave,/spline)
               var[i,j,*] = interpol(cube.var[i,j,*],waveold,wave,/spline)
               dq[i,j,*] = interpol(cube.dq[i,j,*],waveold,wave)
            endfor
         endfor
      ENDIF
      IF cube.ndim eq 2 then begin
         for i=0,cube.ncols-1 do begin
            dat[*,i] = interpol(cube.dat[*,i],waveold,wave,/spline)
            var[*,i] = interpol(cube.var[*,i],waveold,wave,/spline)
            dq[*,i] = interpol(cube.dq[*,i],waveold,wave)
         endfor
      ENDIF
      IF cube.ndim eq 1 then begin
         dat = interpol(cube.dat,waveold,wave,/spline)
         var = interpol(cube.var,waveold,wave,/spline)
         dq = interpol(cube.dq,waveold,wave)
      endif
      ibd = where(dq gt 0.01 AND dq ne 1,ctbd)
      if ctbd gt 0 then begin
         print,'IFSF_VELCOR: Interpolating DQ; ',string(ctbd,format='(I0)'),$
            ' values > 0.01 set to 1.'
         dq[ibd] = 1
      endif

   endif else begin
      print,'IFSF_VELCOR: Must select either VACCOR or HELIOCOR.'
      stop
   endelse

   ; new wavleength solution
   wds = string(cube.wavedim,format='(I0)')
   if datext gt 0 then begin
      newheader_phu = header.phu
      sxaddpar,newheader_phu,'CRVAL'+wds,crval
      sxaddpar,newheader_phu,'CRPIX'+wds,crpix
      sxaddpar,newheader_phu,'CDELT'+wds,cdelt
   endif
   newheader_dat = header.dat
   newheader_var = header.var
   newheader_dq = header.dq
   sxaddpar,newheader_dat,'CRVAL'+wds,crval
   sxaddpar,newheader_dat,'CRPIX'+wds,crpix
   sxaddpar,newheader_dat,'CDELT'+wds,cdelt
   sxaddpar,newheader_var,'CRVAL'+wds,crval
   sxaddpar,newheader_var,'CRPIX'+wds,crpix
   sxaddpar,newheader_var,'CDELT'+wds,cdelt
   if dqext ne -1 then begin
      sxaddpar,newheader_dq,'CRVAL'+wds,crval
      sxaddpar,newheader_dq,'CRPIX'+wds,crpix
      sxaddpar,newheader_dq,'CDELT'+wds,cdelt
   endif
 
   ;  Write output
   if datext gt 0 then begin
      writefits,outfile,cube.phu,newheader_phu
      writefits,outfile,dat,newheader_dat,/append
   endif else begin
      writefits,outfile,dat,newheader_dat
   endelse
   writefits,outfile,var,newheader_var,/append
   if dqext ne -1 then $
      writefits,outfile,dq,newheader_dq,/append

end
