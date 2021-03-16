import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np


degree_li   = np.array([  5.99,   206.11,  111.78,   14.34,    14.45,   111.53,   81.66,    5.89,   15.01,   208.99]) # R.A.の値を入力(単位はdeg)
distance_li = np.array([358.776, 1467.72, 9132.48, 6457.968, 6457.968, 9132.48, 5544.72, 9132.48, 6457.968, 1467.72]) # 距離を入力(単位はly)

fig = plt.figure(figsize=(10, 10))
fig.add_subplot(111).set_aspect('equal', adjustable='box')
fig.patch.set_alpha(0)
# 画像の読み込み
img = mpimg.imread("ssc2008-10b1_Lrg.jpg")
plt.imshow(img)
# for degree, distance in zip(degree_li,distance_li):
plt.scatter(distance_li/45.6*np.cos(np.radians(-degree_li-90))+1497, distance_li/45.6*np.sin(np.radians(-degree_li-90))+2074,s=1,c="r")

plt.grid(0)
plt.axis("off")

print("SAVE")
plt.savefig("銀河系.jpg",dpi=300,bbox_inches='tight', pad_inches=0)
import subprocess
subprocess.run("銀河系.jpg",shell=True)