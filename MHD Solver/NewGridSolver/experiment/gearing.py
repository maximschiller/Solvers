# Calculating the gear size required for a reasonable potentiometer reading based on the estimated thrust
# calculate the required torque to overcome rotational inertia
# calculate the arc based on gear size to determine how much the pot will rotate in degrees
# currently going with a 1:1 gear ratio --> could try 2:1 and 4:1
# NOTE: all references to radius are diameter

import numpy as np
import matplotlib.pyplot as plt

#################################################################

# main parameters - need to check these values
hThruster = 0.4 # height of the thruster moment arm
FThruster = 10*10**-6 # Estimated force produced by thruster
rMaingear = 4/100 # in meters
rPotgear = 0.1*rMaingear

## calculate torque required
# require that the torque produced by the thruster overcomes the moment of inertia of the potentiometer
Tmainbar = hThruster*FThruster

# moment of inertia calculation of potentiometer
mpot = 0.011 # kg
rpot = 0.5/100
alpha = 10 # angular acceleration
inertiapot = 0.5*mpot*rpot**2

Tinertia = inertiapot*alpha

print("Torque produced by thruster: " + str(Tmainbar))
print("Torque required to overcome because of inertia: " + str(Tinertia))

points = 50
x = np.linspace(0,0.1,points)
yinertia = np.ones_like(x)*Tinertia
ymainbar = np.zeros_like(x)
for i in range(0,points):
    # print(x[i])
    ymainbar[i] = (Tmainbar/hThruster)*x[i]

# plt.plot(x, yinertia)
# plt.plot(x, ymainbar)
# plt.show()


#################################################################

## calculate angle that the potentiometer will rotate
# use the relation that the arc length rotated is equal for both gears
# assume the angle that the thruster rotates is 5 degrees
# r1 is the mainbargear and r2 is the pot gear

r1 = rMaingear
r2 = rPotgear
for x in range(0,300,1):
    theta2 = (r1*x)/r2
    if theta2 >= 300:
        print("Max angle of rotation of main bar: " + str(x) + " degrees")
        break


#################################################################
# Assume that the torque of gears can overcome the potentiometer inertia
# arc length is simply defined as s = r*theta
# therfore for two gears who travel the same arc length --> r1*theta1 = r2*theta2
# given that we want the angle the potentiometer can travel to be 300 degrees we can determine the radius ratio of the gears

points = 50
theta2 = 300
r2r1 = np.linspace(0.1,1,50) # the ratio of r2/r1
theta1 = np.zeros_like(r2r1) # initialising theta2


r1 = 6 ## assume a value for radius of main bar in cm
r2 = np.zeros_like(r2r1)
for i in range (0,50):
    theta1[i] = r2r1[i]*theta2
    r2[i] = (r1*theta1[i])/theta2

# plt.plot(r1r2, theta2/theta1, 'ko')
# plt.plot(r1r2, r2, 'bo')
# plt.minorticks_on()
# plt.grid()
# plt.show()

# plt.plot(theta1, r2, 'ko')
# plt.minorticks_on()
# plt.grid()
# plt.show()


# realistically only going to get 5 degrees rotation from main bar
# vary theta2 to see what sort of r2 values i can get
points = 50
angles = 30

theta1 = np.linspace(1,300,angles)
theta2 = 100 # this is the angle of rotatin id ideally like to achieve

r1 = np.linspace(10,100,points) # radius of the main bar gear
r2 = np.zeros((angles, points))

# calculating the radius of the pot gear for a range of main bar gears at every angle
for j in range(0,angles):
    for i in range(0,points):
        r2[j][i] = (r1[i]*theta1[j])/theta2
        

# for j in range(0,angles):
#     plt.plot(r1,r2[j][:])

# # plt.plot(r1[i], r2[j][i])
# plt.minorticks_on()
# plt.grid()
# plt.show()


# creating a 3d contour plot
theta2 = 300
theta1 = np.linspace(1,300,points)

X, Y = np.meshgrid(theta1, r1)
Z = (X*Y)/theta2



# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.minorticks_on
# ax.contour3D(X, Y, Z, 50, cmap='binary')
# ax.set_xlabel('theta1 (deg)')
# ax.set_ylabel('r1 (mm)')
# ax.set_zlabel('r2 (mm)');
# # plt.show()



# calculating r2 after analysing contour plot
# also plotting the resultant torque, recall that T = Fr --> T1/r1 = T2/r2 and have already assumed main bar T
r1 = 37.5 # radius of main bar gear, guess since this is already printed
theta1 = np.linspace(1,300, points)
theta2 = 300 # this is ideally achieved for max sensitivity
T1 = Tmainbar

r2 = np.zeros_like(theta1)
T2 = np.zeros_like(r2)
for i in range(0,points):
    r2[i] = (r1*theta1[i])/theta2

    T2[i] = (T1/r1)*(r1*theta1[i])/theta2


# fig, ax1 = plt.subplots(figsize=(8, 8))
# ax1.plot(theta1, r2,'k')
# ax2 = ax1.twinx()
# ax2.plot(theta1, T2)
plt.plot(theta1,r2)
plt.minorticks_on()
plt.grid()
# plt.show()



#################################################################
# methodology used in solidworks
# 1. assume that the main bar gear will only rotate by a certain amount
# 2. assume a radius for the pot gear that will overcome the inertial torque
# 3. use the ration r1*theta1 = r2*theta2 and use the fact that the required rotation for the pot gear is 300 degrees, i.e. max sensistivity for rotation

# estimated radius of rotation - note that these are diameters
thetaest = 40/300 # want 300 degrees of rotation in pot
d2 = 35 # desired diameter of pot gear

d1 = d2/thetaest
print("Diameter of main bar gear: " + str(d1))

# calculate the required number of teeth for main bar gear
module = 1.5
nteeth = d1/module
print("Number of teeth required: " + str(nteeth))




#################################################################
# GEAR RATIOS USED
# gear ratio 1: 40/360, potgear radius = 36
# gear ratio 2: 80/360, potgear radius = 36
# gear ratio 3: 90/360, potgear radius = 36
# gear ratio 4: 100/360, potgear radius = 42