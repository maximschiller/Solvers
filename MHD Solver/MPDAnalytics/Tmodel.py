import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
import math
import constants

ra = np.arange(0.01, 0.5, 0.01) # radius of the cathode in m
rc = np.arange(0.01, 0.5, 0.01) # radius of the anode in m
ra_center = 0.5 * (ra[:-1] + ra[1:])
rc_center = 0.5 * (rc[:-1] + rc[1:])

RA, RC = np.meshgrid(ra_center, rc_center)

mu0 = constants.magperm()
# Thrust Maecker model
Tsf_Jsq = (mu0/(4*math.pi))*(np.log(RA/RC) + 3/4)
  
# surface plot for Tsf/J^2
fig = plt.figure()

# Uses a 3d prjection as model is supposed to be 3D
axes = fig.add_subplot(projection='3d')

# Plots the three dimensional data consisting of x, y and z 
axes.plot_surface(RA, RC, Tsf_Jsq) 

axes.set_title('3D Parametric Plot')

# Set axes label
axes.set_xlabel('ra', labelpad=5)
axes.set_ylabel('rc', labelpad=5)
axes.set_zlabel('T/Jsq', labelpad=5)

# show command is used to visualize data plot   
plt.show() 


# pcolormesh needs the pixel edges for x and y
# and with default flat shading, Z needs to be evaluated at the pixel center
plot = plt.pcolormesh(ra, rc, Tsf_Jsq, cmap='RdBu', shading='flat')

# contour needs the centers
cset = plt.contour(RA, RC, Tsf_Jsq, cmap='gray')
plt.clabel(cset, inline=True)

plt.colorbar(plot)
plt.show()