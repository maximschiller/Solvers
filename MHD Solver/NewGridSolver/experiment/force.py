## Script to calculate how to adjust parameters in order to achieve desired torque
# given that the force balance equation is --> Fthrustcalc = (Wthruster + Wbar - Wcw*(hcw/hthruster))*np.sin(theta1calc) 
# Useful for when experimenting


# Constant variables
Wthruster = 1
Wbar = 1
hthruster = 1


# Variables that can be changed
Fthrust = 1
Wcw = 1
hcw = 1
theta1 = 1 # note that this is the max angle that can be achieved