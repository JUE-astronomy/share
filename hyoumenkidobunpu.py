from PIL import Image, ImageFilter
import matplotlib.pyplot as plt
import numpy as np
import glob

x_list = []
limli = []
y_list = []
for imagefile in glob.glob(".\\daen\\*")[:10]:
    print(imagefile)
    im2 = Image.open(imagefile)
    # print(im1.getpixel((0, 0)))
    # print(im2.getpixel((0, 0)))
    # print(im1.format, im1.size, im1.mode,im1.getextrema(),im1.getpixel((0, 95)))
    # print(im2.format, im2.size, im2.mode,im2.getextrema(),im2.getpixel((0, 0)))
    c_x,c_y,b_color = [],[],[]
    x_max,y_max,v_max,v_sum = 0,0,0,0
    for x in range(96):
        for y in range(96):
            v = im2.getpixel((x, y))
            v_sum += v
            if v > v_max:
                v_max = v
                x_max = x
                y_max = y

    np.arctan2(y_max,x_max),np.sqrt(x_max**2+y_max**2)
    for x in range(96):
        for y in range(96):
            c_x.append(x-x_max)
            c_y.append(y-y_max)
            b_color.append(im2.getpixel((x, y)))


    b_color = np.array(b_color)
    ax = plt.subplot(1,2,1,polar=True)
    radii = np.sqrt(np.array(c_x)**2 + np.array(c_y)**2)
    theta = np.arctan2(np.array(c_y),np.array(c_x))
    plt.scatter(theta,radii,c=b_color,cmap="Blues")
    ax.grid(0)
    
    sli = []
    limit = 0
    for r in np.linspace(0,80,100):
        s = []
        for n in range(len(radii)):
            if abs(radii[n]) < r:
                s.append(b_color[n])
        if sum(s*0) < v_sum/2:
            sli.append(np.mean(s))
            x_list.append(r)
            y_list.append(np.mean(s))
            limit = r
    limli.append(limit)
    plt.title(str(limit))
    ax = plt.subplot(1,2,2)
    # plt.scatter(range(len(sli)),sli,s=0.05)
    plt.yscale('log')
    plt.ylim(0,300)
    # plt.savefig("subaru_color.png")
    # plt.savefig("subaru_color.png")
    # plt.pause(0.1)
    Image.open(imagefile).close()
    plt.close()
plt.subplot(2,1,1)
plt.scatter(np.array(x_list)**1,y_list)
plt.yscale("log")
plt.subplot(2,1,2)
plt.scatter(np.array(x_list)**(1/2),y_list)
plt.yscale("log")
plt.show()
for name in limli:
    print(name)