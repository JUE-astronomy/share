############# liblaryのimport ################
import warnings

import astropy.io.fits
import numpy as np

warnings.simplefilter('ignore')

# 基本的には////の中を変更する。
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
## 読み込むfileの指定
gas_temp_li   = ['L134N_NH3_Temp.fits', "Orion-B_NH3_Temp.fits", 'N7538_Fuzui_NH3_Temp.fits', 'M17SWex_NH3_Temp.fits', 'M170do_NH3_Temp.fits', 'N7538_NH3_Temp.fits', 'DR21_NH3_Temp.fits', 'W28_NH3_Temp.fits', 'M17SW_NH3_Temp.fits', 'Orion-KL_NH3_Temp.fits']
for nh3 in range(len(gas_temp_li)):
    gas_temp_li[nh3] = "..\\NH3_Temp\\" +(".").join(gas_temp_li[nh3].split(".")[:-1])+"_header_pls.fits"

## 保存する名前を設定
filenames     = ["L134N", "Orion-B", "NGC7538_Fuzui", "M17SWex", "M17_0do", "NGC7538", "DR21", "W28_Fuzui", "M17SW", "Orion-KL"]
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

for gas_temp,filename in zip(gas_temp_li,filenames):
    hdu   = astropy.io.fits.open(gas_temp)[0]
    t_rot = hdu.data
    t_kin = t_rot/(1-t_rot/41.5*np.log10(1+1.1*np.exp(15.7/t_rot)))

    hdu2 = astropy.io.fits.PrimaryHDU(data=t_kin)
    # 元のfitsにあったヘッダーを入れる
    for h in list(hdu.header.keys()):
        try:
            hdu2.header[h] = hdu.header[h]
            print(h,hdu2.header[h])
        except:
            print(h,"Error")

    hdu2.writeto(filename+".fits", overwrite=True)
print("\n","fin.")
