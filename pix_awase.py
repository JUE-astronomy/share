import os
import glob
import shutil
from astropy.io import fits
import numpy as np


#### input
print("INPUT")
obj         = raw_input("object_name (ex. NGC7538) \n   : ")
regrid_list = raw_input('fitsfiles(XXX.fits,YYY.fits...) \n OR directory \n   : ').replace("\\","/").split(',') # XXX.fits,YYY.fits

if len(regrid_list) == 1 and os.path.isdir(regrid_list[0]):
    regrid_list =  glob.glob(regrid_list[0]+"/*.fits")

print('===================================================')
result_list = []
naxis_list  = []
for regrid in regrid_list:
    hdulist = fits.open(regrid)
    try:
        hdu  = hdulist[0]
        x = hdu.header['NAXIS1']
        y = hdu.header['NAXIS2']
    except:
        hdu1 = hdulist[1]
        x = hdu1.header['NAXIS1']
        y = hdu1.header['NAXIS2']
    naxis_list.append(x+y)

min_val   = np.inf
for index,val in enumerate(naxis_list):
    if val < min_val:
        min_val   = val
        min_index = index
template = regrid_list[min_index]

for regrid in regrid_list:
    file  = regrid.split("/")[-1]
    directory = regrid[:-len(file)]
    image = '.image'
    pre   = 'regrid_'
    li = []
    hdulist = fits.open(regrid)
    hdu  = hdulist[0]
    data = hdu.data
    header1 = hdu.header
    try:
        x = hdu.header['NAXIS1']
        y = hdu.header['NAXIS2']
    except:
        hdu1 = hdulist[1]
        x = hdu1.header['NAXIS1']
        y = hdu1.header['NAXIS2']
        data = hdu1.data
        header1 = hdu.header+hdu1.header
    try:
        wave = hdu.header['WAVELEN']
    except:
        wave = hdu.header['WAVELNTH']
    hdu.header['OBJECT'] = obj

    ### saturate delete
    for i in range(y):
        for j in range(x):
            v = data[i][j]
            if v == np.nan :
                v = np.nan
            elif v <= 0:
                v = np.nan
            li.append(v)
    data = np.reshape(li,[y,x]) # reshpe(x*y)
    head = fits.PrimaryHDU(data = data)
    head.header = header1
    filename = directory+obj+"_"+str(wave)+".fits"
    head.writeto(filename, overwrite=True)
    print("Saturated delete has finished.")

    ### regrid
    result = directory+pre+obj+"_"+str(wave)+".fits"
    ### CASAtasks
    importfits(fitsimage=filename, imagename=filename + image)
    importfits(fitsimage=template, imagename=template + image)
    imregrid(imagename=filename + image, output= result+image,template=template + image)
    exportfits(imagename=result+image, fitsimage= result, overwrite=True)
    result_list.append(result)
    print("regrid has finished.")

### create new folder
os.mkdir(obj+"_match")
for name in result_list:
    shutil.move(name,obj+"_match")
print("FINISHED!")