# File which calculates thrust, ISP and efficiency from arduino text file output
# calculate thrust using the rotation of the pot which determines the moment arm
# calculate ISP and efficiency using the accelerometer

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os


# main parameters
g = 9.81
Wthruster = 1
Wbar = 1
Wcw = 1
hcw = 1
hthruster = 1

# NOTE: need to use radius for this
r1 = 1
r2 = 1

# test parameters
mdot = 1
voltage = 1
vexhaust = 1
capacitance = 1

#################################################################

## IMPORT VALUES FROM DATA FILES
# loop through each logger file
# NOTE: structure of files are csv --> millis, angle, accelx, accely, accelz
# specify the directory path
directory = r"C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\solver\NewGridSolver\experiment\data\test"

# create an empty list to store the data
logs = []

# loop over all the log files from log00 to log99
for i in range(100):
    filename = os.path.join(directory, f"LOGGER{i:02d}.csv")  # format the filename with leading zeros and join with directory path
    try:
        # read in the CSV file and append the data to the list
        data = pd.read_csv(filename, header=None)
        logs.append(data.values)
    except FileNotFoundError:
        pass  # ignore missing files

# print the data from the third column of the first file
# print(logs[0][:,2])

len_logs = len(logs)

# storing the values of milliseconds to a list
milli = []
millinoheader = []
for i, data in enumerate(logs):
    millidat = data[:,0]
    milli.append(millidat)

for millidat in milli:
    millinoheaderdat = millidat[1:].astype(float)
    millinoheader.append(millinoheaderdat)

print(millinoheader)
#################################################################

# calculating the force
# force is simply the force body analysis and the max angle achieved
# note that the force is a ratio due to the gear sizing
# need to account for the actual angle travelled by thruster using the arc length relation
# r1*theta1 = r2*theta2

# looping through all the log files and getting the data
# create a new array to store the 3rd column values from each log
theta2_list = []

# loop over the data in logs and assign the 2nd column values to the new array
for i, data in enumerate(logs):
    theta2dat = data[:,1]
    theta2_list.append(theta2dat)

# print(theta2_list)

# # Example of how to convert data for operation
# # Create an empty list to store the float dtype arrays
# float_arrays = []

# # Iterate through each object dtype array in angles_list
# for theta2dat in theta2_list:
#     # Convert the object dtype array to a float dtype array and exclude the first element ("angle")
#     float_array = theta2dat[1:].astype(float)
#     # Append the float dtype array to the list
#     float_arrays.append(float_array)

# # Convert the list of float dtype arrays to a 2D NumPy array
# float_array_2d = np.array(float_arrays)

# print(float_array_2d)

theta2 = []
theta1 = []
Fthrust = []
angles_list = []
# perform a mathematical operation on each column of each log file
for theta2dat in theta2_list:
    # value read by potentiometer - this is theta2
    angles = theta2dat[1:].astype(float)
    angles_list.append(angles)
    degrees = np.deg2rad(angles)   
    theta1calc = (r2/r1)*degrees
    Fthrustcalc = (Wthruster + Wbar - Wcw*(hcw/hthruster))*np.sin(theta1calc)  


    # append the calculated values to their respective lists
    theta2.append(degrees)
    theta1.append(theta1calc)
    Fthrust.append(Fthrustcalc)

# Printing outputs for checking
# print(angles_list)
print(Fthrust)

# # plot the data for each log file
fig = plt.figure()

for i in range(1, len(millinoheader)):
    plt.plot(millinoheader[i], Fthrust[i], label=f"log{i:02d}")

# add labels and legend
plt.xlabel("milliseconds")
plt.ylabel("Thrust")
plt.legend()
plt.show()

#################################################################

# plotting the varying acceleration

# storing the values of acceleration to a list
accelxheader = []
accelyheader = []
accelzheader = []
for i, data in enumerate(logs):
    accelxdat = data[:,2]
    accelxheader.append(accelxdat)
    accelydat = data[:,2]
    accelyheader.append(accelydat)
    accelzdat = data[:,2]
    accelzheader.append(accelzdat)

# removing headers and converting data to a numpy array type
accelx = []
accely = []
accelz = []

for accelxdat in accelxheader:
    accelxnoheaddata = accelxdat[1:].astype(float)
    accelx.append(accelxnoheaddata)
for accelydat in accelyheader:
    accelynoheaddata = accelydat[1:].astype(float)
    accely.append(accelynoheaddata)
for accelzdat in accelzheader:
    accelznoheaddata = accelzdat[1:].astype(float)
    accelz.append(accelznoheaddata)

# create three subplots
fig, axs = plt.subplots(3, 1, figsize=(8, 10))

# loop over the logs and plot acceleration values against theta1 on each subplot
for i in range(1, len(millinoheader)):
# for i, data in enumerate(logs):
    # milli = data.iloc[:, 0].values
    # accelx = data.iloc[:, 2].values
    # accely = data.iloc[:, 3].values
    # accelz = data.iloc[:, 4].values
    
    # axs[0].plot(milli, accelx, label=f"log {i}")
    # axs[1].plot(milli, accely, label=f"log {i}")
    # axs[2].plot(milli, accelz, label=f"log {i}")

    axs[0].plot(millinoheader[i], accelx[i], label=f"log {i}")
    axs[1].plot(millinoheader[i], accely[i], label=f"log {i}")
    axs[2].plot(millinoheader[i], accelz[i], label=f"log {i}")


# add axis labels and legend
axs[0].set_ylabel("acceleration x")
axs[1].set_ylabel("acceleration y")
axs[2].set_ylabel("acceleration z")
axs[2].set_xlabel("milliseconds")
axs[0].legend()
axs[1].legend()
axs[2].legend()

# show the plot
plt.show()



#################################################################

# calculating ISP/efficiency/impulse bit/thrust-weight ratio
# use/dependent on values of propellant mass flow rate, exhaust velocity and input power (voltage)
# use papers to determine how to calculate efficiency