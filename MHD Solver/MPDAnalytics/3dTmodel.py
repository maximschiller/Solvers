import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import constants
import math

# Creating dataset
ra,rc,J = np.meshgrid(np.arange(0.01, 0.5, 0.01), 
                    np.arange(0.01, 0.5, 0.01), 
                    np.arange(0, 10000, 500))

mu0 = constants.magperm()
# Thrust Maecker model
Tsf = (mu0/(4*math.pi))*(np.log(ra/rc) + 3/4)*J**2

# Creating plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
sc = ax.scatter3D(ra, rc, J, alpha = 0.8, c = Tsf, cmap = plt.get_cmap('jet'), marker = '.')
plt.colorbar(sc)
# Set axes label
ax.set_xlabel('ra', labelpad=5)
ax.set_ylabel('rc', labelpad=5)
ax.set_zlabel('T/Jsq', labelpad=5)
plt.show()