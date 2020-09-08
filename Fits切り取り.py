import astropy.io.fits
from astropy.nddata import Cutout2D


infile = input("ファイル名(.fits) : ").split(",") # 読み込むファイル名

for k in range(0,len(infile)):
    print("ファイルトリミング "+str(k+1)+" 個目")

    hdulist = astropy.io.fits.open(infile[k])
    hdu = hdulist[0]
    data = hdu.data
    header1 = hdu.header
    x = hdu.header['NAXIS1']
    y = hdu.header['NAXIS2']
    print("・")

    Cx = 505.69408   # 中心の座標
    Cy = 565.99588
    position = (Cx,Cy)
    Sx = 499.99305   # サイズ
    Sy = 766.65601
    size = (Sx,Sy)     # pixels
    cutout = Cutout2D(data, position, size)
    hdu = astropy.io.fits.PrimaryHDU(data = cutout.data)
    print("・")
    hdu.header = header1
    hdu.header['NAXIS1'] = Sx
    hdu.header['NAXIS2'] = Sy

    print("・")
    
    cutout_filename = "sample"+str(k)+".fits" # 出力ファイル名
    hdu.writeto(cutout_filename, overwrite=True)
    print("ファイルトリミング "+str(k+1)+" 個目終了")
    print(" ")
print("Finish")