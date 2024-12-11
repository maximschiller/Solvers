# Script which calculates the required power supply/power bank size for testing - i.e. number and size of capacitors
# determine electrical/thermal conductivities/resistance of materials

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint


## calculate required current for a desired thrust force
mu = (10**-7)*4*np.pi
F = 100*10**-3
ra = 10/1000
rc = 5/1000

I = np.sqrt((F*4*np.pi)/(mu*(np.log(ra/rc + 0.75))))
print("Current requirent in amps: " + str(I))


#####################################################
## power bank parameters
ncapsseries = 5
ncapsparallel = 5
numcaps = ncapsseries*ncapsparallel
# numcaps = 2
V = 1350
Vcap = 400
IndivCcap = 10*10**-6

Ccap = ncapsparallel*(1/(ncapsseries/(IndivCcap)))
Rres = 1/(ncapsparallel/(ncapsseries*0.1))
Lind = 1/(ncapsparallel/(ncapsseries*(10*10**-9)))
# Ccap = (0.9*10**-6)*2
# Rres = 1/(2/0.1)
# Lind = 1/(2/(10*10**-9))

print('Capacitance: ' + str(Ccap))
print('Resistance = ' + str(Rres))
print('Inductance = ' + str(Lind))

RLsq = (Rres/(2*Lind))**2
LC = 1/(Lind*Ccap)
print('(R/L)^2 = ' + "{:e}".format(RLsq) + ' LC = ' + "{:e}".format(LC)) # need to ensure RLsq is greater than LC

## Plot the discharge current vs time
def f(u,x):
    return (u[1], (1/Lind)*(V - (1/Ccap)*u[0] - Rres*u[1]))

time = np.linspace(0,10*10**-6,100)

charge, current = odeint(f, [0, 0], time).T
plt.plot(time, current, 'blue', linewidth = 2)
plt.xlabel('time (s)')
plt.ylabel('current (A)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.grid(True)
plt.savefig(r'C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs' + '/dischargeA.png')
plt.show()

#####################################################
## Estimated thrust from current level

Imax = np.max(current)
print('Peak current achieved (kA) = ' + str(Imax/1000))
Fmax = ((mu*Imax**2)/(4*np.pi))*(np.log(ra/rc) + 0.75)
print('Maximum theoretical thrust (N) = ' + str(Fmax))

Ftheo = []
for i in current:
    F = ((mu*i**2)/(4*np.pi))*(np.log(ra/rc) + 0.75)
    Ftheo.append(F)
    
plt.plot(time, Ftheo, 'blue', linewidth = 2)
plt.xlabel('time (s)')
plt.ylabel('Force (N)')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.grid(True)
plt.savefig(r'C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\Thesis B figs' + '/pulseforce.png')
plt.show()


#####################################################
## Calculate energy storage and power

# E = numcaps*0.5*IndivCcap*Vcap**2
E = Ccap*1350**2
print('Energy storage in caps = ' + str(E))
sec = 2*10**-6 # operational timeof capacitor discharge
W = E/sec
print('Wattage of capacitors (kW) = ' + str(W/1000))


#####################################################
## Calculating bleeder resistor value

Rbleeder = 5000000

#time constant and discharge time
tau = Rbleeder*IndivCcap
tdischarge = 10*tau
print('time taken for capacitor discharge (s) = ' + str(tdischarge))

#Power dissipated in bleeder resistor
Pcap = Vcap**2/Rbleeder
Ptot = numcaps*Pcap
print('Power dissipated in each cap (W) = ' + str(Pcap))
print('Power dissipated in whole circuit (W) = ' + str(Ptot))


#Skin depth effect - not necessary for perf board