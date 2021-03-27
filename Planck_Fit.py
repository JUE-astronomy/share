import glob
import math
import os
import shutil
import warnings

import astropy.io.fits
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.constants import *
from scipy.optimize import curve_fit

warnings.simplefilter('ignore')


def planck(x, A, T):
    x = 10**(np.array(x)-6)
    return np.log10(A*(2*h*c**2/x**5*1/(math.e**(h*c/(k*x*T))-1)))


def ConvUnit(cdelt1, cdelt2):
    return 1/(abs(cdelt1)*abs(cdelt2)*((2*np.pi)/360)**2)*1e-6


n = 0
no_sati_number = 5  # FITに必要な点の数を指定

xli = np.array(np.linspace(0, 4, 10000))
bairitu_initial  = 10
temp_initial     = 7
bairitu_mintoMax = [1e-4,1e+5]
temp_mintoMax    = [5.8, 100]

object_dir_list = input("ディレクトリ名(複数ある場合は\",\"で区切る) : ").split(",")
dir_path = os.path.dirname(os.path.abspath(__file__))+"/"
for object_dir in object_dir_list:
    fits_list,infiles        = [], []
    naxis1_list, naxis2_list = [], []
    temperature_list    = []
    r_squared_list           = []
    for f in glob.glob(dir_path+object_dir+"/*.fits"):
        print("Load :", f)
        infiles.append(f)
    if infiles == []:
        raise ValueError("\""+dir_path+"\"内に目的のディレクトリがありません")
    try:
        shutil.rmtree(dir_path+"Fit_"+object_dir)
    except:
        pass
    os.mkdir(dir_path+"Fit_"+object_dir)
    for file in infiles:
        fits = astropy.io.fits.open(file)
        naxis1_list.append(fits[0].header["NAXIS1"])
        naxis2_list.append(fits[0].header["NAXIS2"])
        fits_list.append(fits[0])
    naxis1 = min(naxis1_list)
    naxis2 = min(naxis2_list)

    for n2 in range(naxis2):
        for n1 in range(naxis1):
            wavelist, vallist = [], []
            for fits in fits_list:
                v = fits.data[n2][n1]
                if not (np.isnan(v) or v == -1*np.inf or v == np.inf) and v >= 0:
                    try:
                        wavelen = fits.header["WAVELEN"]
                    except:
                        wavelen = fits.header["WAVELNTH"]
                    if wavelen == 70 or wavelen == 100 or wavelen == 160:
                        cdelt1, cdelt2 = fits.header.get(
                            'CDELT1'), fits.header.get('CDELT2'),
                        v = v*ConvUnit(cdelt1, cdelt2)
                    elif wavelen == 250:
                        v = v*91
                    elif wavelen == 350:
                        v = v*51.4
                    elif wavelen == 500:
                        v = v*23.6
                    wavelist.append(wavelen)
                    vallist.append(v)
            if ((no_sati_number <= len(wavelist)) and (no_sati_number <= len(vallist))):
                wavelist = np.log10(np.array(wavelist))
                vallist = np.log10(np.array(vallist))
                param_bounds = ((bairitu_mintoMax[0], temp_mintoMax[0]), (bairitu_mintoMax[1], temp_mintoMax[1]))
                syokichi = [bairitu_initial,temp_initial]
                param, cov = curve_fit(planck, wavelist, vallist, p0=syokichi, maxfev=100000, bounds=param_bounds)   # p0
                
                # R^2の値を求める
                residuals = vallist - planck(wavelist, param[0], param[1])
                rss = np.sum(residuals**2)  # residual sum of squares = rss
                tss = np.sum((vallist-np.mean(vallist))**2)
                r_squared = 1 - (rss / tss)
                temperature_list.append(param[1])
                r_squared_list.append(r_squared)
                """
                ############### プロットする ################
                plt.close()
                plt.plot(xli, planck(xli, param[0], param[1]),ls="--",label="Dust_Temp = "+str(round(param[1],3))+" K\n"+"R$^2$ = "+str(round(r_squared,3)))
                plt.scatter(wavelist, vallist)
                plt.legend()
                plt.xlim(0, 3)
                plt.ylim(-2, 4)
                plt.pause(0.5)
                """
                
            else:
                temperature_list.append(np.nan)
                r_squared_list.append(np.nan)
            n += 1
            print("\r"+"..."+str('{:.2f}'.format(((n)/(naxis1*naxis2*len(object_dir_list)))*100))+"%",end="")

    # R^2値のexcelファイル作成
    reshaped_r_squared = reversed(np.reshape(r_squared_list, [naxis2, naxis1]))
    df = pd.DataFrame(reshaped_r_squared)
    df.to_excel(dir_path+'Fit_'+object_dir+'/r.xlsx',
                encoding='utf-8', index=False, header=False)

    # dustの温度マップ作製
    reshaped_temperature = reversed(np.reshape(temperature_list, [naxis2, naxis1]))
    df = pd.DataFrame(reshaped_temperature)
    df.to_excel(dir_path+'Fit_'+object_dir+'/Dust_Temperature.xlsx',
                encoding='utf-8', index=False, header=False)
    fits = astropy.io.fits.open(repr(infiles[0])[1:-1])
    data = fits[0].data
    hdu = astropy.io.fits.PrimaryHDU(data=df[::-1])
    hdu.header = fits[0].header
    # 入れたいヘッダ
    hdu.header['BITPIX'] = -64
    hdu.header['NAXIS'] = 2
    hdu.header['NAXIS1'] = len(df.columns)
    hdu.header['NAXIS2'] = len(df)
    hdu.header['BUNIT'] = 'K'
    # 作ったfitsの保存
    hdu.writeto(dir_path+'Fit_'+object_dir+'/Dust_Temperature.fits', overwrite=True)
    
    # fitsを閉じる
    for file in infiles:
        astropy.io.fits.open(file).close()
    del hdu