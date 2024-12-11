#Code for communicating with Arduino and rotary motion sensor
#Sensor is on A0, data sent to Python from Arduino via serial port
#Note we used pip install PySerial to install package
import serial #Package for serial communications
import time #Used for delays and to assign time codes to data readings
import numpy as np #Used for creating vectors, etc. 
import matplotlib.pyplot as plt #For plotting
import csv #For writing data to a csv file
import math
arduino = serial.Serial('COM14', 115200, timeout=.1) #Open connection to Arduino
samples = 300 #We will take this many readings
angle_data = np.zeros(samples) #Creates a vector for our angle data
time_data = np.zeros(samples) #Creates a vector of same length for time
i = 1;
calibrate = 0 #Value to zero potentiometer reading when pendulum is motionless, read from Arduino
while i!=samples:
    data = arduino.readline()[0:-2].decode('utf-8')
    if data:
        angle_data[i] = float(data)
        time_data[i] = time.perf_counter()
        print(angle_data[i])
        i = i + 1
        angle_data = (angle_data - calibrate)*math.pi/180
plt.plot(time_data,angle_data,'ro')
plt.axis([0,time_data[samples-1],-math.pi/4,math.pi/4])
plt.xlabel("Time (seconds)")
plt.ylabel("Angle (Radians)")
plt.title("Pendulum Motion - Measured with Arduino and potentiometer")
plt.show()
arduino.close()	    
