"""Does ...

Args:
    files (string):
    xpeaks (int): 
    ypeaks (int): 

Returns:
    None

Note:
"""

from astropy.io import fits

# # check args
# narg=len(sys.argv)

# # should be three (including routine name)
# if narg != 3:
#     print("Usage: python kcwi_masksky_ds9.py <imagename> <regionname>")
#     print("imagename : used for array dimensions and filename purposes, ")
#     print("            must be an _intf image.")
#     print("regionname: name of region file containing ds9 mask regions")
#     print("            (typically a .reg)")
#     exit()

# # read arg values
# imfname=sys.argv[1]
# regfname=sys.argv[2]

files = ['cubea.fits','cubeb.fits','cubec.fits','cubed.fits']
xpeaks = [40.2,38.9,20.7,19.3]
ypeaks = [34.1,34.3,34.5,34.7]
for i in range(len(files)):
    cube = fits.open(files[i], mode='update')
    cube[0].header['XPEAK'] = xpeaks[i]
    cube[0].header['YPEAK'] = ypeaks[i]
    cube.close(output_verify='ignore')
