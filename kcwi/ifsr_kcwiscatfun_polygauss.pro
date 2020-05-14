; docformat = 'rst'
;
;+
;
; Model for fitting scattered light to KCWI 2D data. Options include single
; polynomial (default), two polynomials (one for first half of columns, one 
; for second half), and a gaussian or moffat profile on top of the polynomials.
;
; :Categories:
;    IFSRED/KCWI
;
; :Returns:
;    Model evaluated at columns X.
;    
; :Params:
;    X: in, required, type=dblarr
;      Column numbers.
;    P: in, required, type=dblarr(N)
;      Model parameters. Polynomial parameters are first, followed by
;      Gaussian or Moffat parameters as specified below.
;      
; :Keywords:
;    gauss: in, optional, type=bin
;      Add a gaussian profile to model. Calls function GAUSS1 packaged with
;      MPFIT, with parameter order (MEAN, SIGMA, AREA).
;    moffat: in, optional, type=bin
;      Add a moffat profile to model. Calls function DRT_MOFFAT from DRTOOLS,
;      with parameter order mimicking MPFIT version:
;        A[0] = peak value
;        A[1] = profile center
;        A[2] = alpha = (FWHM/2)/sqrt[2^(1/beta)-1]
;        A[3] = beta = "Moffat index"
;    splitpoly: in, optional, type=bin
;      If set, split column data in half and fit separate polymnomial to each.
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
;    Change History::
;      2018aug14, DSNR, created
;      2020may14, DSNR, documented
;
; :Copyright:
;    Copyright (C) 2018--2020 David S. N. Rupke
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
function ifsr_kcwiscatfun_polygauss, X, P, gauss=gauss, moffat=moffat, $
                                     splitpoly=splitpoly
   if keyword_set(moffat) then begin
      deg = n_elements(p)-5
      ymod = DRT_MOFFAT(X, P[deg+1:deg+4])
   endif else if keyword_set(gauss) then begin
      deg = n_elements(p)-4
      ymod = GAUSS1(X, P[deg+1:deg+3])
   endif else begin
      deg = n_elements(p)-1
      ymod = dblarr(n_elements(X))
   endelse
   if keyword_set(splitpoly) then begin
      deg = (deg-1)/2
      imid = floor(double(n_elements(x))/2d)
      ymod[0:imid-1] += poly(x[0:imid-1],p[0:deg])
      ymod[imid:n_elements(x)-1] += $
         poly(x[imid:n_elements(x)-1],p[deg+1:deg+deg+1])
   endif else begin
      ymod += POLY(X,P[0:deg])
   endelse

   RETURN,ymod
END
