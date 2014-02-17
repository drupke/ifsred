; docformat = 'rst'
;
;+
;
; Plot color-coded maps of a quantity across the GMOS FOV in the native GMOS
; spaxel space (i.e., a grid of hexagons).
;
; :Categories:
;    GMOS
;
; :Returns:
;    Postscript map.
;
; :Params:
;    nslits: in, required, type=integer, default=1
;      Number of GMOS "slits" (1 or 2).
;    z: in, required, type=dblarr(N)
;      Values to map (one per spaxel).
;    x: in, required, type=dblarr(N)
;      IFS x coordinates of spaxels.
;    y: in, required, type=dblarr(N)
;      IFS y coordinates of spaxels.
;    outfile: in, required, type=string
;      Name of postscript file to contain output (not incl. ".eps" extension).
;    
; :Keywords:
;    zran: in, optional, type=dblarr(2)
;      Range for map colors
;    ncbdiv: in, optional, type=integer
;      Number of colorbar divisions
;    cbform: in, optional, type=string
;      Formatting string for colorbar labels
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
;      2010jun01, DSNR, created
;      2014jan31, DSNR, added detailed documentation
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
pro gmos_maps_fiber,nslits,z,x,y,outfile,zran=zran,ncbdiv=ncbdiv,cbform=cbform

  if ~keyword_set(zran) then zran=[min(z),max(z)]
  if ~keyword_set(ncbdiv) then ncbdiv=4
  if ~keyword_set(cbform) then cbform='(D0.2)'

  set_plot,'ps',/copy,/interpolate
  if nslits eq 2 then $
     device,/inches,xsize=9.5,ysize=5.5,xoffset=0.5,yoffset=0.5 $
  else device,/inches,xsize=5.5,ysize=5.5,xoffset=0.5,yoffset=0.5
  device,/encapsulated,/color,bits_per_pixel=8
  !P.thick=2
  !P.charsize=2
  !P.charthick=2

  xybuf=0.02d
  dxran = max(x)-min(x)
  dyran = max(y)-min(y)
  xran = [min(x)-xybuf*dxran,max(x)+xybuf*dxran]
  yran = [min(y)-xybuf*dyran,max(y)+xybuf*dyran]
  dzran = zran[1]-zran[0]

  device,filename=outfile+'.eps'

  cbpos=[0.1,0.09,0.9,0.14]
  cgplot,[0],xran=xran,yran=yran,/xsty,/ysty,/iso,position=[0.1,0.25,0.9,0.95]
  cgloadct,13,/reverse
  for i=0,n_elements(z)-1 do begin
     color = byte((z[i]-zran[0])/(zran[1]-zran[0])*255d)
     usersymbol,'HEXAGON2',size=1.75,_extra={color:color,fill:1}
     oplot,[x[i]],[y[i]],psym=8,color=color
  endfor

  ticknames = string(dindgen(ncbdiv+1)*dzran/ncbdiv - $
                     (dzran - zran[1]),format=cbform)
  cgcolorbar,position=cbpos,divisions=ncbdiv,$
             ticknames=ticknames
  
  device,/close_file

end
