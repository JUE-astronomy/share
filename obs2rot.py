import math
import os
import shutil
import time
import warnings

import astropy.io.fits
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import newton

warnings.simplefilter('ignore')

start = int(time.time())

##### 基本的には////の中を変更する。
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
fitsfile_11 = "C:\\Users\\yamah\\Desktop\\JUEN\\seminar\\Study\\extra_analyze\\auto_derivate_Temperature\\NH3_SerpensSouth\\SerpensSouth_NH3_11.cube.fits"
fitsfile_22 = "C:\\Users\\yamah\\Desktop\\JUEN\\seminar\\Study\\extra_analyze\\auto_derivate_Temperature\\NH3_SerpensSouth\\SerpensSouth_NH3_22.cube.fits"

##### fitsをオープン
hdu_11 = astropy.io.fits.open(fitsfile_11)[0]
hdu_22 = astropy.io.fits.open(fitsfile_22)[0]

##### 輝線の範囲を入力
nh311_main_width =  5000, 10000
nh311_sub1_width = 13000, 18000
nh311_sub2_width = -2000,  3000
nh322_main_width =  5000, 10000

##### 3sigma
yuui = 3

#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

savedir_11 = os.path.dirname(fitsfile_11)+"\\emission_11"
savedir_22 = os.path.dirname(fitsfile_11)+"\\emission_22"

for save in [savedir_11, savedir_22]:
    try:
        shutil.rmtree(save)
    except FileNotFoundError:
        pass
    os.mkdir(save)

frequency_11_li = np.arange(hdu_11.header["CRVAL3"]-hdu_11.header["CRPIX3"]*hdu_11.header["CDELT3"],hdu_11.header["CRVAL3"]+hdu_11.header["NAXIS3"]*hdu_11.header["CDELT3"]-hdu_11.header["CRPIX3"]*hdu_11.header["CDELT3"],hdu_11.header["CDELT3"])
frequency_22_li = np.arange(hdu_22.header["CRVAL3"]-hdu_22.header["CRPIX3"]*hdu_22.header["CDELT3"],hdu_22.header["CRVAL3"]+hdu_22.header["NAXIS3"]*hdu_22.header["CDELT3"]-hdu_22.header["CRPIX3"]*hdu_22.header["CDELT3"],hdu_22.header["CDELT3"])

xmin_11,xmax_11 = np.nanmin(frequency_11_li), np.nanmax(frequency_11_li)
xmin_22,xmax_22 = np.nanmin(frequency_22_li), np.nanmax(frequency_22_li)

number, nh311_main_li, nh311_sub_li, nh322_main_li = 0,[],[],[]
m1_li, s1_li, s2_li, m2_li = [],[],[],[]

for frq in range(len(frequency_11_li)):
    if   nh311_sub1_width[0] <= frequency_11_li[frq] <= nh311_sub1_width[1]:
        s1_li.append(frq)
    elif nh311_main_width[0] <= frequency_11_li[frq] <= nh311_main_width[1]:
        m1_li.append(frq)
    elif nh311_sub2_width[0] <= frequency_11_li[frq] <= nh311_sub2_width[1]:
        s2_li.append(frq)
for frq in range(len(frequency_22_li)):
    if   nh322_main_width[0] <= frequency_22_li[frq] <= nh322_main_width[1]:
        m2_li.append(frq)

for x in range(hdu_11.header["NAXIS1"]):
    for y in range(hdu_11.header["NAXIS2"]):
        ##### 各点(x,y)のスペクトルを取る
        d_11 = hdu_11.data[:,y,x]
        d_22 = hdu_22.data[:,y,x]

        ##### スペクトルを一次関数で近似。
        a_11,b_11 = np.polyfit(frequency_11_li, d_11, 1)
        a_22,b_22 = np.polyfit(frequency_22_li, d_22, 1)

        ##### データから近似直線を引き算(BASE LINEのつもり)。
        d_11 = d_11 - (a_11*frequency_11_li+b_11)
        d_22 = d_22 - (a_22*frequency_22_li+b_22)

        if np.nanmax(d_11[m1_li[0]:m1_li[-1]]) > max([np.nanmax(d_11[s1_li[0]:s1_li[-1]]), np.nanmax(d_11[s2_li[0]:s2_li[-1]])]) and np.nanmax(d_11[m1_li[0]:m1_li[-1]]) > yuui*(np.nanmean(d_11)+np.nanstd(d_11)) and np.mean([np.nanmax(d_11[s1_li[0]:s1_li[-1]]),np.nanmax(d_11[s2_li[0]:s2_li[-1]])]) > yuui*(np.nanmean(d_11) + np.nanstd(d_11)) and np.nanmax(d_22[m2_li[0]:m2_li[-1]]) > yuui*(np.nanmean(d_22) + np.nanstd(d_22)):
            ##### スペクトル図を作成 and 保存
            plt.gca().xaxis.get_major_formatter().set_useOffset(False)
            plt.title("x = "+str(x)+", y = "+str(y))
            plt.plot(frequency_11_li,d_11,lw=.5)
            plt.hlines(yuui*(np.nanmean(d_11) + np.nanstd(d_11)), xmin_11, xmax_11, "tab:orange", linestyles='dashed')
            minimum, np.nanmaximum = np.nanmin(d_11), np.nanmax(d_11)
            for line in [nh311_main_width, nh311_sub1_width, nh311_sub2_width]:
                plt.axes().add_patch(plt.Rectangle([line[0],minimum],line[1]-line[0],np.nanmaximum-minimum,ec="m",ls="dashed",fill=False))
            plt.savefig(savedir_11+"\\x = "+str(x)+", y = "+str(y)+".png")
            plt.close()

            plt.gca().xaxis.get_major_formatter().set_useOffset(False)
            plt.title("x = "+str(x)+", y = "+str(y))
            plt.plot(frequency_22_li,d_22,lw=.5)
            plt.hlines(yuui*(np.nanmean(d_22) + np.nanstd(d_22)), xmin_22, xmax_22, "tab:orange", linestyles='dashed')
            minimum, np.nanmaximum = np.nanmin(d_22), np.nanmax(d_22)
            plt.axes().add_patch(patches.Rectangle([nh322_main_width[0],minimum],nh322_main_width[1]-nh322_main_width[0],np.nanmaximum-minimum,ec="m",ls="dashed",fill=False))
            plt.savefig(savedir_22+"\\x = "+str(x)+", y = "+str(y)+".png")
            plt.close()
            ##### 有意な輝線をlistに追加
            nh311_main_li.append(np.nanmax(d_11[m1_li[0]:m1_li[-1]]))
            nh311_sub_li.append(np.nanmean([np.nanmax(d_11[s1_li[0]:s1_li[-1]]),np.nanmax(d_11[s2_li[0]:s2_li[-1]])]))
            nh322_main_li.append(np.nanmax(d_22[m2_li[0]:m2_li[-1]])) # 上限値を求める場合は nh322_main_li.append(yuui*(np.nanmean(d_11) + np.nanstd(d_11)))
        else:
            ##### 有意でない輝線をnanとしてlistに追加
            nh311_main_li.append(np.nan)
            nh311_sub_li.append(np.nan)
            nh322_main_li.append(np.nan)
        number += 1
        print("\r","進行度 :",round(number/((hdu_11.header["NAXIS1"]*hdu_11.header["NAXIS2"])*3)*100,3),"%   ",end="")

e = math.e
li_ms, tau = [], []
li_ms = np.array(nh311_main_li)/np.array(nh311_sub_li)

for henka in li_ms:
    try:
        tau.append(newton(lambda x: henka*(1-e**(-0.28*x))-(1-e**(-x)),2))
    except:
        tau.append(np.nan)
    number += 1
    print("\r","進行度 :",round(number/((hdu_11.header["NAXIS1"]*hdu_11.header["NAXIS2"])**3)*100,3),"%",end="\r")

log = math.log
Temperature = [] # empty list

for t,m2,m1 in zip(tau,nh322_main_li,nh311_main_li):
    try:
        Temperature.append(-41.5/log(-0.282/t*log(1-m2/m1*(1-e**(-t)))))
    except:
        Temperature.append(np.nan)
    number += 1
    print("\r","進行度 :",round(number/((hdu_11.header["NAXIS1"]*hdu_11.header["NAXIS2"])*3)*100,3),"%",end="\r")

new_hdu = astropy.io.fits.PrimaryHDU(data=np.reshape(Temperature, (hdu_11.header["NAXIS1"],hdu_11.header["NAXIS2"])).T.tolist())
# 元のfitsにあったヘッダーを入れる
for h in list(hdu_11.header.keys()):
    try:
        new_hdu.header[h] = hdu_11.header[h]
        print(h,new_hdu.header[h])
    except:
        print(h,"None")

new_hdu.writeto(os.path.dirname(fitsfile_11)+"\\Gas_Temperature.fits", overwrite=True)

print("\n",(int(time.time())-start)/60,"分")
print("fin.")