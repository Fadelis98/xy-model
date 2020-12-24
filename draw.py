import csv
import numpy as np
from matplotlib import pyplot as plt
from numpy.core.fromnumeric import size


def Bindershow(lowtemp=0.4, hightemp=0.54, num=15):
    binder = []
    s = np.linspace(5, 30, 6)
    tem = list(np.linspace(lowtemp, hightemp, num+1))
    tem.pop()
    for i in s:
        name = "Binder" + str(int(i)) + ".csv"
        with open(name, newline='') as csvfile:
            csvreader = csv.reader(csvfile)
            b = []
            for row in csvreader:
                b.append(abs(float(row[0])))
            binder.append(b)
        plt.ylim((0.6, 1))
        plt.plot(tem, b)
        plt.xlabel("tem")
        plt.ylabel("U")
    plt.savefig("相变点")
    plt.show()


def kcshow(length):
    m = []
    l = []
    t = []
    name = "MLt.csv"
    with open(name, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            m.append(float(row[0]))
            l.append(float(row[1]))
            t.append(float(row[2]))
    m = np.array(m)
    l = np.array(l)
    t = np.array(t)
    y = m*np.float_power(l, 1/8-2)
    x = t*l
    for i in range(len(m)//length):
        plt.scatter(x[length*i:length*(i+1)],
                    y[length*i:length*(i+1)], marker=".")
    plt.xlabel("tL")
    plt.ylabel("ML^{1/8}")
    plt.savefig("标度指数")
    plt.show()


def Magread(lowtemp=0, hightemp=3, num=31):
    MPMag = []
    WFMag = []
    tem = list(np.linspace(lowtemp, hightemp, num))
    for i in tem:
        name1 = "MPmag" + "{:.6f}".format(i) + ".csv"
        name2 = "Wolffmag" + "{:.6f}".format(i) + ".csv"
        with open(name1, newline='') as csvfile:
            csvreader = csv.reader(csvfile)
            m = []
            for row in csvreader:
                m.append(float(row[0]))
            MPMag.append(m)
        with open(name2, newline='') as csvfile:
            csvreader = csv.reader(csvfile)
            m = []
            for row in csvreader:
                m.append(float(row[0]))
            WFMag.append(m)
    for i in range(len(MPMag)):
        plt.figure(figsize=(8,6))
        y1 = MPMag[i]
        y2 = WFMag[i]
        plt.plot(y1, label="MP")
        plt.plot(y2, label="Wolff")
        plt.ylim((0,400))
        plt.legend()
        plt.title(f"T={tem[i]}")
        plt.show()

# 关联函数
def relread(scale):
    name = f"Rel{scale}.csv"
    with open(name, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        head = next(csvreader)
        dnum = len(head)-1
        relOft = []
        relOfd = [[] for i in range(dnum)]
        tlist = []
        for row in csvreader:
            tlist.append(float(row[0]))
            relOft.append(list(map(float, row[1:])))
            for i in range(dnum):
                relOfd[i].append(float(row[i+1]))
    return (relOft, relOfd, tlist)


def drawRelt(relOft, tlist):
    plt.figure(figsize=(16, 9))
    for i in range(len(relOft)):
        y = relOft[i]
        plt.plot(y, label=f"T={tlist[i]}")
        plt.legend()
    plt.show()


def drawReld(relOfd, tlist):
    plt.figure(figsize=(16, 9))
    for i in range(len(relOfd)):
        y = relOfd[i]
        plt.plot(tlist, y, label=f"d={i}")
        plt.legend()
    plt.show()


def readLat(filename):
    lat = []
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            lat.append(list(map(float, row[:-1])))
    return lat


def drawlat(lat, title):
    plt.figure(figsize=(8, 8))
    plt.imshow(lat, cmap="hsv", vmin=0, vmax=2)
    plt.colorbar()
    plt.title(title)
    plt.show()


def readVort(filename):
    data = []
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            data.append(list(map(float, row[:-1])))
    numOfVort2 = len(data[0])
    Vort = [(data[0][i], data[0][i+1]) for i in range(0, numOfVort2, 2)]
    numOfantiVort2 = len(data[1])
    antiVort = [(data[1][i], data[1][i+1])
                for i in range(0, numOfantiVort2, 2)]
    return Vort, antiVort


def drawlatWithVort(lat, title, Vort, antiVort):
    plt.figure(figsize=(8, 8))
    plt.imshow(lat, cmap="hsv", vmin=0, vmax=2)
    plt.colorbar()
    if(len(Vort) != 0):
        vortx, vorty = zip(*Vort)
        plt.scatter(vortx, vorty, marker="x", s=150, c="#000000")
    if(len(antiVort) != 0):
        avortx, avorty = zip(*antiVort)
        plt.scatter(avortx, avorty, marker="*", s=150, c="#FFFFFF")
    plt.title(title)
    plt.show()


def noteall(Tlist,vList,scale):
    for t in Tlist:
        lat = readLat("LatOf{0:.6f}s{1}.csv".format(t,scale))
        plt.figure(figsize=(8, 8))
        plt.imshow(lat, cmap="hsv", vmin=0, vmax=2)
        plt.colorbar()
        for v in vList:
            try:
                vort, antivort = readVort("VortOf{0:.6f}s{1}v{2}.csv".format(t, scale,v))
                if(len(vort) != 0):
                    vortx, vorty = zip(*vort)
                    plt.scatter(vortx, vorty, marker="*", s=v**2*30, c="#000000")
                if(len(antivort) != 0):
                    avortx, avorty = zip(*antivort)
                    plt.scatter(avortx, avorty, marker="*", s=v**2*30, c="#FFFFFF")
            except IOError:
                pass
        plt.title(f"T = {t}")
        plt.show()


