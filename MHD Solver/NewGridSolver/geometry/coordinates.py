import numpy as np
import matplotlib.pyplot as plt
from itertools import zip_longest

# Main script which generates domain text file based on the given function input
# i.e. whether it be condi, diverging, cylindrical or conical

# reading the coordinates into one large file
def read_coordinates(*args, filename):
    with open(filename, 'w') as file:
        for coordinates in zip_longest(*args, fillvalue=None):
            file.write('\t'.join(str(coord) for coord in coordinates) + '\n')

# seperating the coordinates into different text files based on boundary
def extract_coordinates(input_filename, output_filename, col1, col2):
    # Read the input file and split the lines into columns
    with open(input_filename, 'r') as file:
        lines = file.readlines()
        columns = [line.strip().split() for line in lines]

    # Extract the desired columns from the list of columns, ignoring None values
    # col1_values = [float(column[col1]) for column in columns if column[col1] is not None]
    # col2_values = [float(column[col2]) for column in columns if column[col2] is not None]
    col1_values = [float(column[col1]) for column in columns if column[col1] is not None and column[col1] != 'None']
    col2_values = [float(column[col2]) for column in columns if column[col2] is not None and column[col2] != 'None']


    # Write the extracted columns to the output file
    with open(output_filename, 'w') as file:
        for x, y in zip(col1_values, col2_values):
            file.write(f"{x}\t{y}\n")

