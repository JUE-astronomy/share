########## Use: make a sub-image using a CASA task
########## Last Modified: 2018/02/04
########## Author: Tsuyoshi ISHIDA
##### dependencies
from __future__ import print_function, division
import os, sys
from astropy.io import fits

from astropy.wcs import WCS
import matplotlib.pyplot as plt
##### version check
version = sys.version_info.major

##### arguments
# argv = sys.argv
# argc = len(argv)
# if argc != 2:
#     raise SyntaxError('The number of arguments is wrong.')

##### filenames
# fitsname = argv[1]

##### reference
print('================================================================================')
print('SDP.81:')
print('    Band 6 & 7 continuum: cx, cy, x, y = 374.70952, 330.81108, 250, 250')
print('SDP.9:')
print('================================================================================')

##### input
if version == 2:
    fitsname     = raw_input('trimmed : ')
    template     = raw_input('template: ')
else:
    fitsname     = input('trimmed : ')
    template     = input('template: ')



hdu_tmp  = fits.open(template)[0]
cx_world,cy_world = WCS(template).all_pix2world((hdu_tmp.header["NAXIS1"])/2,(hdu_tmp.header["NAXIS2"])/2,0)
w_world,h_world   = abs(hdu_tmp.header["NAXIS1"]*hdu_tmp.header["CDELT1"]),abs(hdu_tmp.header["NAXIS2"]*hdu_tmp.header["CDELT2"])

hdu_fits = fits.open(fitsname)[0]
cx,cy = WCS(fitsname).all_world2pix(cx_world,cy_world,0)
x,y   = w_world/abs(hdu_fits.header["CDELT1"]), h_world/abs(hdu_fits.header["CDELT2"])


fitsname0, fitsext = os.path.splitext(fitsname)
subsfx  = '_sub'    # suffix for subtracted images
imext   = '.image'  # extension for CASA images
outfile = fitsname0 + subsfx + fitsext

bx, by  = (cx - 1) - x / 2, (cy - 1) - y / 2
tx, ty  = (cx - 1) + x / 2 - 1, (cy - 1) + y / 2 - 1

##### CASA tasks
importfits(fitsimage=fitsname, imagename=fitsname0+imext)
imsubimage(imagename=fitsname0+imext, outfile=fitsname0+subsfx+imext,
           box='{}, {}, {}, {}'.format(bx, by, tx, ty), dropdeg=False)
exportfits(imagename=fitsname0+subsfx+imext, fitsimage=outfile, overwrite=True)

##### change header information
data, header = fits.getdata(outfile, header=True)
try:
    del header['ORIGIN']
except KeyError:
    pass

##### write to FITS
try:
    fits.writeto(outfile, data, header, overwrite=True)
except TypeError:
    fits.writeto(outfile, data, header, clobber=True)

##### delete intermediate files
os.system('rm -rf {} {} {} {}'.format(fitsname0+imext, fitsname0+subsfx+imext, '*.last', '*.log'))
