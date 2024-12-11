## calculate values using analytical formula
import numpy as np

## calculate required current for a desired thrust force
mu = (10**-7)*4*np.pi


# Maecker Model
ra_rc = np.linspace(1,2,10)
current = np.linspace(0,20000,500)

T = []
for i in ra_rc:
    for j in current:
        Tcalc = (mu/(4*np.pi))*(np.log(ra_rc) + 0.75)*current**2
        Tcalc.append(T)

        plt.plot(current,T)

plt.show()