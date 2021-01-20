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


def planck(x, a, b, T1, T2):
    return np.log10(planck1(x, a, T1)+planck1(x, b, T2))


def planck1(x, A, T):
    x = 10**(np.array(x)-6)
    return A*(2*h*c**2/x**5*1/(math.e**(h*c/(k*x*T))-1))


def ConvUnit(cdelt1, cdelt2):
    return 1/((abs(cdelt1)*abs(cdelt2)*3600)/3600*((2*np.pi)/360)**2)*1e-6




n = 0
no_sati_number = 10  # FITに必要な点の数を指定


xli = np.array(np.linspace(0, 4, 10000))
bairitu1_initial  = 1
bairitu2_initial  = 1
temp1_initial     = 110
temp2_initial     = 7
bairitu1_mintoMax = [-500, 1000]
bairitu2_mintoMax = [-500, 1000]
temp1_mintoMax    = [30, 1000]
temp2_mintoMax    = [5.8, 100]

object_dir_list = input("ディレクトリ名(複数ある場合は\",\"で区切る) : ").split(",")
dir_path = os.path.dirname(os.path.abspath(__file__))+"/"
for object_dir in object_dir_list:
    fits_list,infiles        = [], []
    naxis1_list, naxis2_list = [], []
    warm_temperature_list    = []
    cold_temperature_list    = []
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
                        v = v*90
                    elif wavelen == 350:
                        v = v*51.4
                    elif wavelen == 500:
                        v = v*23.6
                    wavelist.append(wavelen)
                    vallist.append(v)
            long = 0
            for w in wavelist:
                if 70 <= w:
                    long += 1
            if ((no_sati_number <= len(wavelist)) and (no_sati_number <= len(vallist))) or (5 <= long):
                wavelist = np.log10(np.array(wavelist))
                vallist = np.log10(np.array(vallist))
                param_bounds = ((bairitu1_mintoMax[0], bairitu2_mintoMax[0], temp1_mintoMax[0], temp2_mintoMax[0]), (
                    bairitu1_mintoMax[1], bairitu2_mintoMax[1], temp1_mintoMax[1], temp2_mintoMax[1]))
                syokichi = [bairitu1_initial, bairitu2_initial,
                            temp1_initial, temp2_initial]
                param, cov = curve_fit(planck, wavelist, vallist, p0=syokichi, maxfev=100000, bounds=param_bounds)   # p0
                """
                ############### プロットする ################
                plt.plot(xli, planck(xli, param[0], param[1], param[2], param[3]))
                plt.plot(xli, np.log10(planck1(xli, param[0], param[2])), ls="--", alpha=0.6,label="Warm = "+str(round(param[2],3)))
                plt.plot(xli, np.log10(planck1(xli, param[1], param[3])), ls="--", alpha=0.6,label="Cold = "+str(round(param[3],3)))
                plt.scatter(wavelist, vallist)
                plt.legend()
                plt.xlim(0, 3)
                plt.ylim(-2, 4)
                plt.show()
                """
                # R^2の値を求める
                residuals = vallist - \
                    planck(wavelist, param[0], param[1], param[2], param[3])
                rss = np.sum(residuals**2)  # residual sum of squares = rss
                # total sum of squares = tss
                tss = np.sum((vallist-np.mean(vallist))**2)
                r_squared = 1 - (rss / tss)
                warm_temperature_list.append(param[2])
                cold_temperature_list.append(param[3])
                r_squared_list.append(r_squared)
            else:
                warm_temperature_list.append(np.nan)
                cold_temperature_list.append(np.nan)
                r_squared_list.append(np.nan)
            n += 1
            print("\r"+"..."+str('{:.2f}'.format(((n)/(naxis1*naxis2*len(object_dir_list)))*100))+"%",end="")

    # R^2値のexelファイル作成
    reshaped_r_squared = reversed(np.reshape(r_squared_list, [naxis2, naxis1]))
    df = pd.DataFrame(reshaped_r_squared)
    df.to_excel(dir_path+'Fit_'+object_dir+'/R_Squared.xlsx',
                encoding='utf-8', index=False, header=False)

    # warm dustの温度マップ作製
    reshaped_warm_temperature = reversed(np.reshape(warm_temperature_list, [naxis2, naxis1]))
    df = pd.DataFrame(reshaped_warm_temperature)
    df.to_excel(dir_path+'Fit_'+object_dir+'/Warm_Dust_Temperature.xlsx',
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
    hdu.writeto(dir_path+'Fit_'+object_dir+'/Warm_Dust_Temperature.fits', overwrite=True)

    # cold dustの温度マップ作製
    reshaped_cold_temperature = reversed(np.reshape(cold_temperature_list, [naxis2, naxis1]))
    df = pd.DataFrame(reshaped_cold_temperature)
    df.to_excel(dir_path+'Fit_'+object_dir+'/Cold_Dust_Temperature.xlsx',
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
    hdu.writeto(dir_path+'Fit_'+object_dir+'/Cold_Dust_Temperature.fits', overwrite=True)
    
    # fitsを閉じる
    for file in infiles:
        astropy.io.fits.open(file).close()
    del hdu