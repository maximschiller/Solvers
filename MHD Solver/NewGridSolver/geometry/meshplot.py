import numpy as np
import matplotlib.pyplot as plt

# importing the values for the grid lines
with open(r"C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\solver\NewGridSolver\meshxvals.txt", "r") as infile:
    arrx = [[float(x) for x in line.split()] for line in infile]
with open(r"C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\solver\NewGridSolver\meshyvals.txt", "r") as infile:
    arry = [[float(y) for y in line.split()] for line in infile]

#importing the values for the grid centres
with open(r"C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\solver\NewGridSolver\centrexvals.txt", "r") as infile:
    xc = [[float(y) for y in line.split()] for line in infile]
with open(r"C:\Users\maxim\OneDrive - The University of Sydney (Students)\THESIS\solver\NewGridSolver\centreyvals.txt", "r") as infile:
    yc = [[float(y) for y in line.split()] for line in infile]


plt.plot(arrx, arry, c='k') # use plot, not scatter
plt.plot(np.transpose(arrx), np.transpose(arry),c='k') # add this here
plt.scatter(xc,yc)
plt.show()