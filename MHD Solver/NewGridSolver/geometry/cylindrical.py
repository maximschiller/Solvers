import numpy as np
import matplotlib.pyplot as plt
from coordinates import read_coordinates
from coordinates import extract_coordinates
import os

spacing = 20

# generating top boundary
radius = 15
cyllen = 60
xcyl = np.linspace(0,cyllen,spacing)
ycyl = np.ones_like(xcyl)*radius

xtop = xcyl
ytop = ycyl
xnegtop = xtop
ynegtop = -ytop


# generating bottom boundary
tipstart = 20-np.sqrt(20**2 - 10**2)
tiprad = 10
xbotcathode = np.linspace(0,tipstart,spacing)
ybotcathode = np.ones_like(xbotcathode)*tiprad
xbotcathodetip = np.linspace(tipstart,20,spacing)
mtip = (0-tiprad)/(20-tipstart)
ybotcathodetip = mtip*xbotcathodetip
c = tiprad - ybotcathodetip[0]
ybotcathodetip = mtip*xbotcathodetip + c
xbotexit = np.linspace(20,xtop[-1],spacing)
ybotexit = np.zeros_like(xbotexit)

xbot = np.concatenate([xbotcathode, xbotcathodetip, xbotexit])
ybot = np.concatenate([ybotcathode, ybotcathodetip, ybotexit])
xnegbot = xbot
ynegbot = -ybot


# generating left boundary
yleft = np.linspace(ybotcathode[0],ytop[0],spacing)
xleft = np.zeros_like(yleft)

yleftneg = -yleft
xleftneg = xleft


# generating right boundary
yright = np.linspace(0,ytop[-1],spacing)
xright = np.ones_like(yright)*xtop[-1]

yrightneg = -yright
xrightneg = xright




# plotting
# positive boundaries
plt.plot(xtop,ytop, 'b')
plt.plot(xbot, ybot, 'r')
plt.plot(xleft,yleft, 'k')
plt.plot(xright,yright,'k')

# negative boundaries
plt.plot(xnegtop,ynegtop, 'b')
plt.plot(xnegbot,ynegbot, 'r')
plt.plot(xleftneg,yleftneg,'k')
plt.plot(xrightneg,yrightneg,'k')

plt.savefig(r'C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs' + '/cylindrical.png')

plt.show()


# outputting points to text file
read_coordinates(xtop, ytop, xleft, yleft, xbot, ybot, xright, yright, filename = 'coordinates.txt')


# calling function which outputs xbot,ybot,xright,yright,xtop,ytop,xleft,yleft
xtopout = []
ytopout = []
xleftout = []
yleftout = []
xbotout = []
ybotout = []
xrightout = []
yrightout = []

# seperating each boundary to a seperate text file
with open('coordinates.txt', 'r') as file:
    for line in file:
        columns = line.split()
        if len(columns) >= 8:
            xtopout.append(columns[0])
            ytopout.append(columns[1])
            xleftout.append(columns[2])
            yleftout.append(columns[3])
            xbotout.append(columns[4])
            ybotout.append(columns[5])
            xrightout.append(columns[6])
            yrightout.append(columns[7])

# writing the coordinates to seperate text files
extract_coordinates('coordinates.txt','topbound.txt',0,1)
extract_coordinates('coordinates.txt','leftbound.txt',2,3)
extract_coordinates('coordinates.txt','botbound.txt',4,5)
extract_coordinates('coordinates.txt','rightbound.txt',6,7)