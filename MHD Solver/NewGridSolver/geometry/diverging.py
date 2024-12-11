import numpy as np
import matplotlib.pyplot as plt
from coordinates import read_coordinates
from coordinates import extract_coordinates
import os


spacing = 20 # spacing between points on nozzle

## Assumed variables
# radius = ybell[0]
radius = 10
cyllen = 30
tipstart = 25-np.sqrt(20**2 - 10**2)
tiprad = 5
Rt = 6 # in mm - just an assumption for now

# Using an 80% bell curve as it is most efficient
# Will need to recalculate this angles for an accurate expantion ratio
# Use eps =70 as this is reasonable for a high altitude leo nozzle, same as the golauncher engines
# thetan = 33
# thetae = 7
# using eps = 8
thetan = 25
thetae = 11

print("Nozzle angle", thetan)
print("Exit angle", thetae)

## Nozzle expansion ratio calculation
# Atmospheric pressure at leo
P_leo = 80000 # in pascals
# Pressure in MPD due to magnetic pressure
P_mpd = 8000 # in pascals

# specific heats constant
gam_air = 1.3

Astar = (((gam_air+1)/2)**(1/(gam_air-1))) * ((P_mpd/P_leo)**(1/gam_air)) * np.sqrt(((gam_air+1)/(gam_air-1)) * (1 - (P_mpd/P_leo)**((gam_air-1)/gam_air)))
eps = 1/Astar
print("Area expansion ratio", eps)
eps = 8


## Calculation of nozzle length
# for an 80% bell length
Re = np.sqrt(eps)*Rt
print("Exit nozzle radius", Re)

Ln = 0.8*(((np.sqrt(eps)-1)*Rt)/np.tan(np.deg2rad(15)))
print("Nozzle length", Ln)


## Throat section calculation
# entrant section
theta = np.linspace(-135, -90, spacing)
xent = 1.5*Rt*np.cos(np.deg2rad(theta))
yent = 1.5*Rt*np.sin(np.deg2rad(theta)) + 1.5*Rt + Rt

# exit section
theta = np.linspace(-90,(thetan-90),spacing)
xex = 0.382*Rt*np.cos(np.deg2rad(theta))
yex = 0.382*Rt*np.sin(np.deg2rad(theta)) + 0.382*Rt + Rt

## Bell section calculation
t = np.linspace(0,1,spacing)

# calculating points Nx and Ny by setting exit section formulas to thetan-90
Nx = 0.382*Rt*np.cos(np.deg2rad(thetan-90))
Ny = 0.382*Rt*np.sin(np.deg2rad(thetan-90)) + 0.382*Rt + Rt

# calculating points Ex and Ey
Ex = Ln
Ey = Re

# calculating points Qx and Qy
m1 = np.tan(np.deg2rad(thetan))
m2 = np.tan(np.deg2rad(thetae))
C1 = Ny - m1*Nx
C2 = Ey - m2*Ex
Qx = (C2-C1)/(m1-m2)
Qy = (m1*C2-m2*C1)/(m1-m2)

# Bell section
xbell = ((1-t)**2)*Nx + 2*(1-t)*t*Qx + (t**2)*Ex
ybell = ((1-t)**2)*Ny + 2*(1-t)*t*Qy + (t**2)*Ey
zbell = np.zeros_like(xbell)


## Just for bell nozzle
diffx = cyllen-xbell[0]
diffy = radius-ybell[0]
# diffx = 0
# diffy = 0
with open(r'solidworksgeom/diverging_cad.txt', 'w') as f:
    # Write header row
    # f.write('x\ty\n')
    
    # Write data rows
    for i in range(len(xbell)):
        f.write('{}\t{}\t{}\n'.format(xbell[i]+diffx, ybell[i]+diffy, zbell[i]))


# generating top boundary
xcyl = np.linspace(0,cyllen,spacing)
ycyl = np.ones_like(xcyl)*radius

xtop = np.concatenate([xcyl, xbell+diffx])
ytop = np.concatenate([ycyl, ybell+diffy])
xnegtop = xtop
ynegtop = -ytop


# generating bottom boundary
xbotcathode = np.linspace(0,tipstart,spacing)
ybotcathode = np.ones_like(xbotcathode)*tiprad
xbotcathodetip = np.linspace(tipstart,20,spacing)
mtip = (0-tiprad)/(20-tipstart)
ybotcathodetip = mtip*xbotcathodetip
c = tiprad - ybotcathodetip[0]
ybotcathodetip = mtip*xbotcathodetip + c
xbotexit = np.linspace(20,xbell[-1]+diffx,spacing)
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


## Plotting
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

plt.savefig(r'C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs' + '/diverging.png')

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