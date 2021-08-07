import os
from plyfile import PlyData, PlyElement
import numpy as np
import pandas as pd
import bats
import matplotlib.pyplot as plt
import time
import scipy.spatial.distance as distance
import csv
# load functions that will record time
from comp_function import time_BATS_rips_flags
from comp_function import time_BATS_updates_rips
from comp_function import time_dionysus_rips
from comp_function import time_gudhi_rips
from comp_function import time_ripser_rips
from comp_function import time_toplayer_rips



def read_ply_to_numpy_arr(file_path):

    plydata = PlyData.read(file_path) # read file
    data = plydata.elements[0].data # read data
    data_pd = pd.DataFrame(data) # Convert to DataFrame, because DataFrame can parse structured data
    data_np = np.zeros(data_pd.shape, dtype=np.float32) # Initialize the array of stored data
    property_names = data[0].dtype.names # read the name of the property
    for i, name in enumerate(property_names): # Read data according to property, so as to ensure that the data read is the same data type.
        data_np[:, i] = data_pd[name]

    return data_np[:,:3] # only use the first 3 columns

# Now ready to run
# write into files
header = ['standard_basis','standard_no_basis','clearing','compression','update',
'topologylayer','gudhi','dionysus','ripser'] # set header of csv file


flags = [
    (bats.standard_reduction_flag(), bats.compute_basis_flag()),
    (bats.standard_reduction_flag(),),
    (bats.standard_reduction_flag(), bats.clearing_flag()),
    (bats.standard_reduction_flag(), bats.compression_flag())
]
labels = [
    "standard w/ basis",
    "standard w/ no basis",
    "standard w/ clearing",
    "standard w/ compression",
]

# csv file name
n_sample = 400 # number of sample points 
filename = 'bunny_record{}_H1.csv'.format(n_sample)

with open(filename, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # read and transform into numpy array
    full_X = read_ply_to_numpy_arr('bunny/data/bun000.ply')
    n_total_rows = full_X.shape[0]

    # set parameter
    dmax_complex = 2 # maximum dimension of complex, will compute PH until one dimension lower
    data_into_file = [] # create a list used to write into file
    n_repeats = 20
    for i in range(n_repeats):

        # Uniform randomly generate 100 points
        random_indices = np.random.choice(n_total_rows, size=n_sample, replace=False)
        X = full_X[random_indices, :]
        X2 = X + X * 0.01 * np.random.randn(X.shape[0],X.shape[1])

        # write row
        record_line = [] # record for each line

        # BATs Method with Flags
        for flag, label in zip(flags, labels):
            record_line.append(time_BATS_rips_flags(X2, dmax_complex, flag))

        # BATs update method
        record_line.append(time_BATS_updates_rips(X, X2, dmax_complex))

        # topologylayer
        record_line.append(0) # topologylayer is slow not run
        # record_line.append(time_toplayer_rips(X2,dmax_complex-1))

        # Gudhi
        record_line.append(time_gudhi_rips(X2, dmax_complex))
        
        # dionysus
        record_line.append(0) # dionysus slow not run
        # record_line.append(time_dionysus_rips(X2, d)) # dionysus slow not run
        
        # ripser
        # Take care, maxdim is Maximum homology dimension
        # dmax_complex is the complex dimension, need to mius one
        record_line.append(time_ripser_rips(X2, maxdim = dmax_complex-1)) 
        
        # add row into data
        data_into_file.append(record_line)

    # write multiple rows
    writer.writerows(data_into_file)
